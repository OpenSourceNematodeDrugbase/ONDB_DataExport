"""
wbpCelegansOrthologues.py (v1)

Purpose:
    For ONE parasite species (species_code), fetch C. elegans orthologue
    information from WormBase ParaSite BioMart and reduce to the
    "best" C. elegans orthologue per parasite Gene stable ID, based on
    highest "% identity", with deterministic tie-breaking.

Contract (mirrors wbpHumanOrthologues.py):

    - header="1" XML
    - central fetch via queryWbpBiomart.fetch_biomart_with_validation
    - single-species input (no commas)
    - Polars internally, Pandas at boundary
    - no renaming of BioMart labels (we treat BioMart header as truth)
    - one row per "Gene stable ID" in the final Pandas output
"""

from __future__ import annotations

import time
from typing import Set

import pandas as pd
import polars as pl

from queryWbpBiomart import fetch_biomart_with_validation  # central validated fetch


# ===== CONFIGURATION (match populateGeneList.py / wbpHumanOrthologues.py) =====

REQUEST_TIMEOUT = 120
RATE_LIMIT_DELAY = 0.5
MAX_RETRIES = 3


# ===== Canonical column labels from header="1" discovery =====
# These MUST match the BioMart headers exactly.

PARASITE_GENE_ID_COL = "Gene stable ID"
PARASITE_GENOME_COL = "Genome project"

CE_GENE_ID_COL = "Caenorhabditis elegans (PRJNA13758) [WS290] gene stable ID"
CE_GENE_NAME_COL = "Caenorhabditis elegans (PRJNA13758) [WS290] gene name"
CE_ORTHO_TYPE_COL = "Homology type"
CE_PCT_ID_COL = "% identity"
CE_PCT_ID_R1_COL = "Caenorhabditis elegans (PRJNA13758) [WS290] % identity"

# Minimal required columns for this module to function
EXPECTED_CELEGANS_ORTHO_COLUMNS: Set[str] = {
    PARASITE_GENE_ID_COL,
    CE_GENE_ID_COL,
    CE_GENE_NAME_COL,
    CE_ORTHO_TYPE_COL,
    CE_PCT_ID_COL,
    CE_PCT_ID_R1_COL,
}


# ===== Contract-compliant C. elegans orthologue XML =====

CELEGANS_ORTHO_XML_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "wbps_gene" interface = "default" >
		<Filter name = "species_id_1010" value = "{species_code}"/>
		<Filter name = "biotype" value = "protein_coding"/>
		<Attribute name = "production_name_1010" />
		<Attribute name = "wbps_gene_id" />
		<Attribute name = "caelegprjna13758_gene" />
		<Attribute name = "caelegprjna13758_gene_name" />
		<Attribute name = "caelegprjna13758_orthology_type" />
		<Attribute name = "caelegprjna13758_homolog_perc_id" />
		<Attribute name = "caelegprjna13758_homolog_perc_id_r1" />
	</Dataset>
</Query>
"""


# ===== Helpers =====

def _validate_single_species(species_code: str) -> str:
    if not isinstance(species_code, str) or not species_code.strip():
        raise ValueError("[error] species_code must be a non-empty string")
    if "," in species_code:
        raise ValueError(f"[error] Exactly one species_code expected (no commas): {species_code!r}")
    return species_code.strip()


def _build_xml(species_code: str) -> str:
    """
    Interpolate species_code into the fixed XML template.
    This is the ONLY modification to the XML besides header="1".
    """
    species_code = _validate_single_species(species_code)
    return CELEGANS_ORTHO_XML_TEMPLATE.format(species_code=species_code)


def _fetch_with_retry(
    xml_query: str,
    expected_columns: Set[str],
    timeout: int = REQUEST_TIMEOUT,
) -> pl.DataFrame:
    """
    Retry wrapper around the central BioMart fetcher, mirroring the pattern
    used in wbpHumanOrthologues.py and populateGeneList.py.
    """
    last_err = None
    for attempt in range(MAX_RETRIES):
        try:
            if attempt > 0:
                time.sleep(RATE_LIMIT_DELAY * (2 ** attempt))
            df = fetch_biomart_with_validation(
                xml_query,
                expected_columns=expected_columns,
                timeout=timeout,
            )

            # Debug output for header troubleshooting
            print(f"[DEBUG] C. elegans orthologues — ACTUAL COLUMNS: {df.columns}")
            print(f"[DEBUG] C. elegans orthologues — EXPECTED COLUMNS: {expected_columns}")

            return df
        except Exception as e:
            last_err = e
            print(f"[DEBUG] C. elegans orthologues — ERROR: {e}")
            if attempt == MAX_RETRIES - 1:
                break
            print(f"Attempt {attempt + 1} failed, retrying...: {e}")

    raise RuntimeError(
        f"C. elegans orthologues BioMart query failed after {MAX_RETRIES} attempts: {last_err}"
    )


# ===== Biology-preserving reduction (Polars internal) =====

def _reduce_to_best_celegans_orthologue(df: pl.DataFrame) -> pl.DataFrame:
    """
    For each parasite Gene stable ID, keep the C. elegans orthologue row
    with the highest % identity, breaking ties deterministically
    by original row order (mirrors pandas groupby().idxmax behaviour).

    Uses CE_PCT_ID_COL ("% identity") as the selection metric, just like the
    human module uses "% identity".
    """
    # Preserve original order to emulate idxmax "first occurrence"
    df = df.with_row_count("row_nr")

    # Coerce % identity to float and fill nulls with 0.0
    df = df.with_columns(
        pl.col(CE_PCT_ID_COL).cast(pl.Float64, strict=False).fill_null(0.0)
    )

    # Sort by (Gene stable ID, % identity desc, original row order asc)
    df_best = (
        df.sort(
            [PARASITE_GENE_ID_COL, CE_PCT_ID_COL, "row_nr"],
            descending=[False, True, False],
        )
        .unique(subset=[PARASITE_GENE_ID_COL], keep="first")
        .drop("row_nr")
    )

    return df_best


def _add_criteria_and_evidence(df_best: pl.DataFrame) -> pl.DataFrame:
    """
    Add C. elegans orthologue flags and evidence columns.

    Pattern:
        - lacks_celegans_orthologue (bool)
        - lacks_celegans_orthologue_evidence (str)
        - best_celegans_orthologue_lt_40pct_identity (bool)
        - best_celegans_orthologue_lt_40pct_identity_evidence (str)

    Threshold is 40% identity, mirroring the human orthologue module.
    """
    # 1) lacks_celegans_orthologue
    df1 = df_best.with_columns(
        pl.col(CE_GENE_ID_COL).is_null().alias("lacks_celegans_orthologue")
    )

    # Evidence for lacks_celegans_orthologue
    lacks_ev = (
        pl.when(pl.col("lacks_celegans_orthologue"))
        .then(pl.lit("No C. elegans orthologue listed in WormBase ParaSite"))
        .otherwise(
            pl.lit("Has C. elegans orthologue(s) in WormBase ParaSite, the most similar is: ")
            + pl.col(CE_GENE_NAME_COL)
        )
        .alias("lacks_celegans_orthologue_evidence")
    )

    # 2) Threshold-based flag for "weak" C. elegans orthologue
    THRESHOLD = 40.0

    ce_lt_thresh = (
        (pl.col(CE_PCT_ID_COL) < THRESHOLD)
        .alias("best_celegans_orthologue_lt_40pct_identity")
    )

    pct_as_str = pl.col(CE_PCT_ID_COL).cast(pl.Utf8)
    ce_name = pl.col(CE_GENE_NAME_COL).fill_null("")

    ce_lt_thresh_ev = (
        pl.when(pl.col("best_celegans_orthologue_lt_40pct_identity"))
        .then(
            pl.lit("Best C. elegans orthologue in WormBase ParaSite has < 40% identity: ")
            + ce_name + pl.lit(" ") + pct_as_str + pl.lit("%")
        )
        .otherwise(
            pl.lit("Best C. elegans orthologue found in WormBase ParaSite has >= 40% identity: ")
            + ce_name + pl.lit(" ") + pct_as_str + pl.lit("%")
        )
        .alias("best_celegans_orthologue_lt_40pct_identity_evidence")
    )

    # Override evidence when there is no C. elegans orthologue at all
    ce_lt_thresh_ev_final = (
        pl.when(pl.col("lacks_celegans_orthologue"))
        .then(pl.lit("No C. elegans orthologue listed in WormBase ParaSite"))
        .otherwise(pl.col("best_celegans_orthologue_lt_40pct_identity_evidence"))
        .alias("best_celegans_orthologue_lt_40pct_identity_evidence")
    )

    out = (
        df1
        .with_columns([lacks_ev, ce_lt_thresh])
        .with_columns([ce_lt_thresh_ev])        # temporary
        .with_columns([ce_lt_thresh_ev_final])  # override when lacks == True
    )

    return out


def _assert_integrity(df: pl.DataFrame) -> None:
    """
    Integrity checks that do not inject new biology:
        - no duplicate 'Gene stable ID' after reduction
    """
    dup_count = df.filter(pl.col(PARASITE_GENE_ID_COL).is_duplicated()).height
    if dup_count > 0:
        raise ValueError(
            f"Duplicate values found in '{PARASITE_GENE_ID_COL}' after C. elegans reduction: {dup_count}"
        )


# ===== Public entry: Polars → Pandas boundary =====

def retrieveCelegansOrthologuesFromWbpBiomart(species_code: str) -> pd.DataFrame:
    """
    Fetch C. elegans orthologue fields (header=1) for ONE parasite species
    using the central BioMart client, then:

        - reduce to one C. elegans orthologue per parasite Gene stable ID
          based on highest % identity
        - add flags/evidence columns
        - assert no duplicate Gene stable ID

    Returns:
        Pandas DataFrame with at least:

            PARASITE_GENE_ID_COL
            CE_GENE_ID_COL
            CE_GENE_NAME_COL
            CE_ORTHO_TYPE_COL
            CE_PCT_ID_COL
            CE_PCT_ID_R1_COL
            lacks_celegans_orthologue
            lacks_celegans_orthologue_evidence
            best_celegans_orthologue_lt_40pct_identity
            best_celegans_orthologue_lt_40pct_identity_evidence
    """
    species_code = _validate_single_species(species_code)
    print(f"[info] Processing species (C. elegans orthologues): {species_code}")

    # 1) Build exact XML and fetch (validated; raises on empty)
    xml_query = _build_xml(species_code)
    df = _fetch_with_retry(xml_query, EXPECTED_CELEGANS_ORTHO_COLUMNS, timeout=REQUEST_TIMEOUT)

    print(f"[info] BioMart (C. elegans orthologues) — received columns: {list(df.columns)}")
    print(f"[info] BioMart (C. elegans orthologues) — received rows: {df.height:,}")

    # 2) Biology-preserving reduction: max % identity per Gene stable ID
    df_best = _reduce_to_best_celegans_orthologue(df)

    # 3) Criteria + evidence columns
    df_out = _add_criteria_and_evidence(df_best)

    # 4) Integrity (no duplicate Gene stable ID)
    _assert_integrity(df_out)

    # 5) Boundary: return Pandas
    return df_out.to_pandas()


# Optional: small CLI for smoke-testing
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="C. elegans orthologues: fetch+reduce (max % identity) + criteria/evidence"
    )
    parser.add_argument("--species-code", required=True)
    args = parser.parse_args()

    out = retrieveCelegansOrthologuesFromWbpBiomart(args.species_code)
    print("OK. Shape:", out.shape)
    print("Columns:", list(out.columns))
