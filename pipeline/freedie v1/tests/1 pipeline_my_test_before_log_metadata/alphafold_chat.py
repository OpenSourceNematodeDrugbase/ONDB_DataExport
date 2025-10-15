#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#this code works - only searchs for genes that has TRUE values of AF

import requests
import pandas as pd
from io import StringIO

MARTSERVICE = "https://parasite.wormbase.org/biomart/martservice"

XML_QUERY = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
  <Dataset name = "wbps_gene" interface = "default" >
    <Filter name = "with_alphafold" excluded = "0"/>
    <Filter name = "species_id_1010" value = "wubancprjeb536"/>
    <Attribute name = "production_name_1010" />
    <Attribute name = "wbps_gene_id" />
  </Dataset>
</Query>
"""

def fetch_wbp_biomart_using_xml(xml_query: str) -> str:
    """POST the XML to BioMart and return TSV text."""
    r = requests.post(MARTSERVICE, data={"query": xml_query}, timeout=60)
    r.raise_for_status()
    return r.text

def run(out_csv="wbancrofti_genes_with_alphafold.csv"):
    # 1) Fetch
    tsv = fetch_wbp_biomart_using_xml(XML_QUERY)

    # 2) Parse (order matches the requested attributes)
    df = pd.read_csv(StringIO(tsv), sep="\t", header=None, names=["production_name", "gene_id"])

    # 3) Clean up: drop duplicates, add flag
    df = df.drop_duplicates().reset_index(drop=True)
    df["alphafold_available"] = True

    # 4) (Optional) add AlphaFold URL via UniProt if you later add UniProt attrs.
    # For now, we only certify presence (True) because we're filtering by with_alphafold.

    # 5) Save + preview
    df.to_csv(out_csv, index=False)
    print(f"Rows: {len(df)}  |  Saved: {out_csv}")
    print(df.head(10).to_string(index=False))

if __name__ == "__main__":
    run()
