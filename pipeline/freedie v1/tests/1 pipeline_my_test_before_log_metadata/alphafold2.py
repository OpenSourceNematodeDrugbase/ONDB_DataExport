#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#This option for W.bancrofti (wubancprjna275548) does return alpha fold structures

import requests
import pandas as pd
from io import StringIO

MARTSERVICE = "https://parasite.wormbase.org/biomart/martservice"

def build_xml(species_key: str, excluded_flag: int) -> str:
    """
    excluded_flag = 0 -> return genes WITH AlphaFold
    excluded_flag = 1 -> return genes WITHOUT AlphaFold
    """
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="parasite_mart" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
  <Dataset name="wbps_gene" interface="default">
    <Filter name="with_alphafold" excluded="{excluded_flag}"/>
    <Filter name="species_id_1010" value="{species_key}"/>
    <Attribute name="production_name_1010" />
    <Attribute name="wbps_gene_id" />
  </Dataset>
</Query>
"""

def fetch_wbp_biomart_using_xml(xml_query: str) -> str:
    r = requests.post(MARTSERVICE, data={"query": xml_query}, timeout=60)
    r.raise_for_status()
    return r.text

def query_genes(species_key: str, excluded_flag: int) -> pd.DataFrame:
    """
    Returns DataFrame with columns: production_name, gene_id
    """
    xml = build_xml(species_key, excluded_flag)
    tsv = fetch_wbp_biomart_using_xml(xml)
    if not tsv.strip():
        # Empty result
        return pd.DataFrame(columns=["production_name", "gene_id"])
    df = pd.read_csv(StringIO(tsv), sep="\t", header=None, names=["production_name", "gene_id"])
    return df.drop_duplicates()

def build_alphafold_table(species_key: str = "wubancprjna275548",
                          out_csv: str = "wbancrofti_alphafold_by_gene_wubancprjna275548.csv") -> pd.DataFrame:
    # A) Genes WITH AlphaFold
    df_with = query_genes(species_key, excluded_flag=0)
    df_with["alphafold_available"] = True

    # B) Genes WITHOUT AlphaFold
    df_without = query_genes(species_key, excluded_flag=1)
    df_without["alphafold_available"] = False

    # C) Union → one row per gene with True/False flag
    # Prefer production_name from whichever side has it (they should match anyway)
    all_df = pd.concat([df_with, df_without], ignore_index=True)
    # If a gene (unexpectedly) appeared in both, True should win
    all_df = (all_df.sort_values("alphafold_available", ascending=False)
                    .drop_duplicates(subset=["gene_id"], keep="first"))

    # D) Pretty order + save
    all_df = all_df[["gene_id", "production_name", "alphafold_available"]].sort_values("gene_id")
    all_df.to_csv(out_csv, index=False)

    # E) Print totals
    total = len(all_df)
    n_true = int(all_df["alphafold_available"].sum())
    n_false = total - n_true
    print(f"Species: {species_key}")
    print(f"Total genes: {total} | With AlphaFold: {n_true} | Without AlphaFold: {n_false}")
    print(f"Saved: {out_csv}")
    print(all_df.head(12).to_string(index=False))
    return all_df

if __name__ == "__main__":
    # Default to the assembly key used in the XML:
    # - Legacy build: "wubancprjeb536"
    # If later we want the newer assembly, swap to e.g. "name of the project"
    build_alphafold_table("wubancprjna275548")
