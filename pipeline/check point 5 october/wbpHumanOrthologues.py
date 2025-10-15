"""
wbpHumanOrthologues.py (Stage 2)
– Biology copied exactly from the original script:
  • For each parasite 'Gene stable ID', keep the row with the highest '% identity'
  • Add four columns:
      - lacks_WBP_human_orthologue (bool)
      - lacks_WBP_human_orthologue_evidence (str)
      - best_WBP_human_orthologue_lt_40pct_identity (bool, threshold = 40)
      - best_WBP_human_orthologue_lt_40pct_identity_evidence (str)
– Engineering dynamics mirror populateGeneList_2.py:
  • header="1" XML
  • central fetch + REQUIRED columns validation (raises on empty)
  • retry/backoff wrapper
  • Polars internal → Pandas at boundary
"""

from __future__ import annotations

import time
from typing import Set
import pandas as pd
import polars as pl

from queryWbpBiomart import fetch_biomart_with_validation  # central validated fetch

# ===== CONFIGURATION (match populateGeneList_2.py) =====
REQUEST_TIMEOUT = 120
RATE_LIMIT_DELAY = 0.5
MAX_RETRIES = 3

# ===== EXPECTED COLUMNS from BioMart (header=1, your source-of-truth) =====
EXPECTED_HUMAN_ORTHO_COLUMNS: Set[str] = {
    "Gene stable ID",
    "Human gene stable ID",
    "Human gene name",
    "Homology type",
    "% identity",
    "Human % identity",
}

# ===== FIXED XML TEMPLATE (exact; header="1") =====
HUMAN_ORTHO_XML_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
  <Dataset name = "wbps_gene" interface = "default" >
    <Filter name = "species_id_1010" value = "{species_code}"/>
    <Filter name = "biotype" value = "protein_coding"/>
    <Attribute name = "wbps_gene_id" />
    <Attribute name = "hsapiens_gene" />
    <Attribute name = "hsapiens_gene_name" />
    <Attribute name = "hsapiens_orthology_type" />
    <Attribute name = "hsapiens_homolog_perc_id" />
    <Attribute name = "hsapiens_homolog_perc_id_r1" />
  </Dataset>
</Query>"""

# ===== Retry wrapper (same style as populateGeneList_2.py) =====
def _fetch_with_retry(xml_query: str, expected_columns: Set[str], timeout: int = REQUEST_TIMEOUT) -> pl.DataFrame:
    last_err = None
    for attempt in range(MAX_RETRIES):
        try:
            if attempt > 0:
                time.sleep(RATE_LIMIT_DELAY * (2 ** attempt))
            df = fetch_biomart_with_validation(xml_query, expected_columns=expected_columns, timeout=timeout)

            # Debug mirrors gene-list module
            print(f"[DEBUG] ACTUAL COLUMNS: {df.columns}")
            print(f"[DEBUG] EXPECTED COLUMNS: {expected_columns}")

            # `fetch_biomart_with_validation` already validates columns & non-empty
            return df
        except Exception as e:
            last_err = e
            print(f"[DEBUG] ERROR: {e}")
            if attempt == MAX_RETRIES - 1:
                break
            print(f"Attempt {attempt + 1} failed, retrying.: {e}")
    raise RuntimeError(f"BioMart query failed after {MAX_RETRIES} attempts: {last_err}")

# ===== Biology-preserving reduction (Polars internal) =====
def _reduce_to_best_orthologue(df: pl.DataFrame) -> pl.DataFrame:
    """
    EXACT behaviour from the original script:
    - fill null '% identity' with 0
    - for each 'Gene stable ID', keep the row with the highest '% identity'
    - tie behaviour mirrors pandas groupby().idxmax(): keep the first occurrence among ties
      (we achieve this by preserving original row order and breaking ties by that order).
    """
    # Preserve original order to emulate pandas idxmax “first occurrence”
    df = df.with_row_count("row_nr")

    # Coerce and fill for selection step
    df = df.with_columns(
        pl.col("% identity").cast(pl.Float64, strict=False).fill_null(0.0)
    )

    # Sort by ('Gene stable ID', '% identity' desc, original order asc) and keep first per gene
    df_best = (
        df.sort(["Gene stable ID", "% identity", "row_nr"], descending=[False, True, False])
          .unique(subset=["Gene stable ID"], keep="first")
          .drop("row_nr")
    )
    return df_best

def _add_criteria_and_evidence(df_best: pl.DataFrame) -> pl.DataFrame:
    """
    EXACT strings and threshold from the original:
      - 40% threshold
      - evidence wording unchanged (“WormBase ParaSite”, “found in WormBase Parasite…”)
    """
    # lacks_WBP_human_orthologue
    df1 = df_best.with_columns(
        pl.col("Human gene stable ID").is_null().alias("lacks_WBP_human_orthologue")
    )

    # Evidence for lacks_WBP_human_orthologue
    lacks_ev = (
        pl.when(pl.col("lacks_WBP_human_orthologue"))
          .then(pl.lit("No human orthologue listed in WormBase ParaSite"))
          .otherwise(
              pl.lit("Has human orthologue(s) in WormBase ParaSite, the most similar is: ") + pl.col("Human gene name")
          )
          .alias("lacks_WBP_human_orthologue_evidence")
    )

    # best_WBP_human_orthologue_lt_40pct_identity
    lt40 = (pl.col("% identity") < 40).alias("best_WBP_human_orthologue_lt_40pct_identity")

    # Evidence for < 40% (use stringifying behaviour like original pandas .astype(str))
    pct_as_str = pl.col("% identity").cast(pl.Utf8)
    gene_name = pl.col("Human gene name").fill_null("")

    lt40_ev = (
        pl.when(pl.col("best_WBP_human_orthologue_lt_40pct_identity"))
          .then(
              pl.lit("Best human orthologue in WormBase ParaSite has < 40% identity: ")
              + gene_name + pl.lit(" ") + pct_as_str + pl.lit("%")
          )
          .otherwise(
              pl.lit("Best human orthologue found in WormBase Parasite has >= 40% identity: ")
              + gene_name + pl.lit(" ") + pct_as_str + pl.lit("%")
          )
          .alias("best_WBP_human_orthologue_lt_40pct_identity_evidence")
    )

    # Override <40 evidence when lacks orthologue (exact original behaviour)
    lt40_ev_final = (
        pl.when(pl.col("lacks_WBP_human_orthologue"))
          .then(pl.lit("No human orthologue listed in WormBase ParaSite"))
          .otherwise(pl.col("best_WBP_human_orthologue_lt_40pct_identity_evidence"))
          .alias("best_WBP_human_orthologue_lt_40pct_identity_evidence")
    )

    out = (
        df1
        .with_columns([lacks_ev, lt40])
        .with_columns([lt40_ev])  # temporary
        .with_columns([lt40_ev_final])  # override when lacks == True
    )
    return out

def _assert_integrity(df: pl.DataFrame) -> None:
    """
    Keep the integrity checks that do NOT inject new biology:
    - No duplicate 'Gene stable ID' after reduction.
    (We drop the old ≥10,000 rows heuristic per our contract.)
    """
    dup_count = df.filter(pl.col("Gene stable ID").is_duplicated()).height
    if dup_count > 0:
        raise ValueError(f"Duplicate values found in 'Gene stable ID' after reduction: {dup_count}")

# ===== Public entry: Polars → Pandas boundary =====
def retrieveHumanOrthologuesFromWbpBiomart(species_code: str) -> pd.DataFrame:
    """
    Fetch human-orthologue fields (header=1) for ONE species using the central client,
    then reproduce the original biology: keep max '% identity' per gene and add criteria/evidence.
    Returns Pandas so run_v01.py can merge to export.csv.
    """
    if not isinstance(species_code, str) or not species_code.strip():
        raise ValueError("species_code must be a non-empty string")
    if "," in species_code:
        raise ValueError("This function processes ONE species at a time.")

    species_code = species_code.strip()
    print(f"Processing species: {species_code}")

    # 1) Build exact XML and fetch (validated; raises on empty like populateGeneList_2.py)
    xml_query = HUMAN_ORTHO_XML_TEMPLATE.format(species_code=species_code)
    df = _fetch_with_retry(xml_query, EXPECTED_HUMAN_ORTHO_COLUMNS, timeout=REQUEST_TIMEOUT)
    print(f"[info] BioMart (human orthologues) — received columns: {list(df.columns)}")
    print(f"[info] BioMart (human orthologues) — received rows: {df.height:,}")



    # 2) Biology-preserving reduction
    df_best = _reduce_to_best_orthologue(df)

    # 3) Criteria + evidence columns (exact wording/threshold)
    df_out = _add_criteria_and_evidence(df_best)

    # 4) Integrity (no dup Gene stable ID). No 10k-size assertion anymore.
    _assert_integrity(df_out)

    # 5) Boundary: return Pandas
    return df_out.to_pandas()

# ----- Optional: small CLI for smoke-testing when the query returns rows -----
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Human orthologues: fetch+reduce (max % identity) + criteria/evidence")
    parser.add_argument("--species-code", required=True)
    args = parser.parse_args()
    out = retrieveHumanOrthologuesFromWbpBiomart(args.species_code)
    print("OK. Shape:", out.shape)
    print("Columns:", list(out.columns))
