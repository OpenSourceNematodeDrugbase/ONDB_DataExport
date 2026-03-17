from queryWbpBiomart import *


# Expected columns match EXACT BioMart header=1 labels
EXPECTED_HUMAN_ORTHO_COLUMNS = {
    "Gene stable ID",
    "Human gene stable ID",
    "Human gene name",
    "Homology type",
    "% identity",
    "Human % identity",
}


def print_biomart_header_labels(species_code: str) -> None:
    xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
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

    # Empty set = no strict validation; just return the DataFrame so we can see the exact labels.
    df = fetch_biomart_with_validation(xml_query, expected_columns=EXPECTED_HUMAN_ORTHO_COLUMNS)


    # df = fetch_wbp_biomart(xml_query, engine="polars", timeout=120)
    print("\n=== BioMart header labels (source of truth) ===")
    for i, c in enumerate(df.columns, 1):
        print(f"{i:2d}. {c}")
    print("\nPaste-ready Python list for EXPECTED_INPUT_COLUMNS:")
    print(repr(df.columns))



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-code", default="wbancrofti_prjeb427")
    args = parser.parse_args()
    print_biomart_header_labels(args.species_code)

#run with python 1.py --species-code wbancrofti_prjeb427

# === BioMart header labels (source of truth) ===
#  1. Gene stable ID
#  2. Human gene stable ID
#  3. Human gene name
#  4. Homology type
#  5. % identity
#  6. Human % identity

# Paste-ready Python list for EXPECTED_INPUT_COLUMNS:
# ['Gene stable ID', 'Human gene stable ID', 'Human gene name', 'Homology type', '% identity', 'Human % identity']