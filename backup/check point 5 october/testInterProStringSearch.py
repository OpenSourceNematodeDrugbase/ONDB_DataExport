# testInterProStringSearch.py — Step 2 (string-search logic added)
# Contract:
# - EXACT XML (header="1"); only interpolate species_code
# - species_code single-species validation
# - Header discovery prints
# - Minimal required-columns validation
# - Polars internally; Pandas at boundary
# - No renaming of BioMart headers; only adds <test_name> and <test_name>_evidence

from __future__ import annotations

import polars as pl
import pandas as pd

# Central fetcher (use whichever one exists in your repo)
try:
    from queryWbpBiomart import fetch_wbp_biomart_using_xml_polars
except ImportError:
    from queeryingwbpbiomart import fetch_wbp_biomart_using_xml_polars


def _validate_single_species(species_code: str) -> str:
    if not isinstance(species_code, str) or not species_code.strip():
        raise ValueError("[error] species_code must be a non-empty string")
    if "," in species_code:
        raise ValueError(f"[error] Exactly one species_code expected (no commas): {species_code!r}")
    return species_code.strip()


def testInterProStringSearch(species_code: str, search_string: str, test_name: str) -> pd.DataFrame:
    """
    Returns a Pandas DataFrame with columns:
      - "Gene stable ID"
      - <test_name> (bool)
      - <test_name>_evidence (str)

    Behaviour:
      - Build "interpro_annotation" = "InterPro ID - InterPro description"
      - Mark a row as match if interpro_annotation regex-contains `search_string`
      - Per gene: if any match -> True and evidence lists matched annotations (unique, sorted, "; " joined)
                  else -> False and a negative evidence sentence
    """
    species_code = _validate_single_species(species_code)

    # === EXACT XML (header="1"); attributes match your original ===
    xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
              
        <Dataset name = "wbps_gene" interface = "default" >
        <Filter name = "species_id_1010" value = "{species_code}"/>
        <Filter name = "biotype" value = "protein_coding"/>
        <Attribute name = "wbps_gene_id" />
        <Attribute name = "interpro_id" />
        <Attribute name = "interpro_description" />
        </Dataset>
    </Query>
    """

    # Fetch (Polars)
    df: pl.DataFrame = fetch_wbp_biomart_using_xml_polars(xml_query)

    # Header discovery
    cols = df.columns
    print(f"[info] Biomart header discovered ({len(cols)} columns):")
    for i, c in enumerate(cols, start=1):
        print(f"  [{i:02d}] {c}")

    # Minimal required set
    required = {"Gene stable ID", "InterPro ID", "InterPro description"}
    missing = [c for c in required if c not in cols]
    if missing:
        raise RuntimeError(f"[error] Missing required columns: {missing}. Check XML and Biomart headers.")

    # If no rows at all, return empty shell with the right schema
    if df.height == 0:
        empty = pl.DataFrame(
            {"Gene stable ID": pl.Series([], dtype=pl.Utf8),
             test_name: pl.Series([], dtype=pl.Boolean),
             f"{test_name}_evidence": pl.Series([], dtype=pl.Utf8)}
        )
        return empty.to_pandas()

    # Build annotation and match flag
    df2 = (
        df
        .with_columns(
            (pl.col("InterPro ID") + pl.lit(" - ") + pl.col("InterPro description")).alias("interpro_annotation")
        )
        .with_columns(
            pl.when(pl.col("interpro_annotation").is_not_null())
              .then(pl.col("interpro_annotation").str.contains(search_string))  # regex
              .otherwise(pl.lit(False))
              .alias("_match")
        )
    )

    # Aggregate to one row per gene
    grouped = (
        df2
        .group_by("Gene stable ID", maintain_order=True)
        .agg([
            pl.col("_match").any().alias(test_name),
            pl.col("interpro_annotation").filter(pl.col("_match")).alias("_hits"),
        ])
        .with_columns([
            pl.when(pl.col(test_name))
              .then(
                  pl.lit("Encodes protein with InterPro domain(s): ")
                  + pl.col("_hits").list.unique().list.sort().list.join("; ")
              )
              # Negative sentence is GPCR-oriented because this module is used for GPCR in your pipeline
              .otherwise(pl.lit("Encodes protein that lacks InterPro GPCR domains"))
              .alias(f"{test_name}_evidence"),
        ])
        .select(["Gene stable ID", test_name, f"{test_name}_evidence"])
    )

    # Integrity: one row per gene
    if grouped.height != grouped.select(pl.col("Gene stable ID")).unique().height:
        raise AssertionError("[error] Duplicate 'Gene stable ID' rows detected in module output.")

    return grouped.to_pandas()
