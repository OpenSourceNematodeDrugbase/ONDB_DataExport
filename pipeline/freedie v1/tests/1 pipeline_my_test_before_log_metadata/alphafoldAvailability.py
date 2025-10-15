
#this code coded on the 16 sept - linked to runpipeline2.py

# alphafoldAvailability.py
import requests
import pandas as pd
from io import StringIO

MARTSERVICE = "https://parasite.wormbase.org/biomart/martservice"

def _build_xml(species_keys: str, excluded_flag: int) -> str:
    """
    excluded_flag = 0 -> return genes WITH AlphaFold
    excluded_flag = 1 -> return genes WITHOUT AlphaFold
    species_keys: single key or comma-separated list (e.g. "wubancprjeb536,wubancprjna275548")
    """
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="parasite_mart" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
  <Dataset name="wbps_gene" interface="default">
    <Filter name="with_alphafold" excluded="{excluded_flag}"/>
    <Filter name="species_id_1010" value="{species_keys}"/>
    <Attribute name="production_name_1010" />
    <Attribute name="wbps_gene_id" />
  </Dataset>
</Query>
"""

def _fetch_biomart(xml_query: str) -> str:
    r = requests.post(MARTSERVICE, data={"query": xml_query}, timeout=60)
    r.raise_for_status()
    return r.text

def _query_block(species_keys: str, excluded_flag: int) -> pd.DataFrame:
    """
    Returns: DataFrame with columns: production_name, Gene stable ID
    """
    xml = _build_xml(species_keys, excluded_flag)
    tsv = _fetch_biomart(xml)
    if not tsv.strip():
        return pd.DataFrame(columns=["production_name", "Gene stable ID"])
    df = pd.read_csv(StringIO(tsv), sep="\t", header=None, names=["production_name", "Gene stable ID"])
    return df.drop_duplicates()

def queryAlphaFoldAvailability(genomes: str) -> pd.DataFrame:
    """
    Main entrypoint.
    Input:
        genomes: WBPS species key(s), comma-separated (e.g., "trtricprjeb535,wubancprjna275548")
    Output:
        DataFrame with columns:
          - Gene stable ID
          - alphafold_available (bool)
          - production_name (to report “how many genomes scanned”)
    """
    df_with = _query_block(genomes, excluded_flag=0)
    df_with["alphafold_available"] = True

    df_without = _query_block(genomes, excluded_flag=1)
    df_without["alphafold_available"] = False

    # Union; if any accidental duplicates slip through, True wins
    all_df = pd.concat([df_with, df_without], ignore_index=True)
    all_df = (all_df.sort_values("alphafold_available", ascending=False)
                    .drop_duplicates(subset=["Gene stable ID"], keep="first"))

    # Nice order
    all_df = all_df[["Gene stable ID", "production_name", "alphafold_available"]].sort_values("Gene stable ID")

    # Console stats (optional—leave printing to pipeline)
    genes_total = all_df["Gene stable ID"].nunique()
    genes_with = (all_df["alphafold_available"] == True).sum()
    genes_without = genes_total - genes_with
    genomes_total = all_df["production_name"].nunique()

    print(f"[AlphaFold] GENES scanned -> With: {genes_with:,} | Without: {genes_without:,} | Total: {genes_total:,}")
    print(f"[AlphaFold] GENOMES scanned (distinct production_name): {genomes_total:,}")

    return all_df
