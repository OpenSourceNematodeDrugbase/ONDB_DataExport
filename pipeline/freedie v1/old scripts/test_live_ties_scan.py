import argparse
import pandas as pd
from queryWbpBiomart import fetch_wbp_biomart_using_xml  # old fetcher (returns as-is, no sort)  # noqa

# This XML matches the old human-orthologue module's attributes.
XML_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
  <Dataset name = "wbps_gene" interface = "default" >
    <Filter name = "species_id_1010" value = "{species}"/>
    <Filter name = "biotype" value = "protein_coding"/>
    <Attribute name = "wbps_gene_id" />
    <Attribute name = "hsapiens_gene" />
    <Attribute name = "hsapiens_gene_name" />
    <Attribute name = "hsapiens_orthology_type" />
    <Attribute name = "hsapiens_homolog_perc_id" />
    <Attribute name = "hsapiens_homolog_perc_id_r1" />
  </Dataset>
</Query>
"""

def main(species: str, sample_n: int):
    # 1) Fetch raw table via the old client (it auto-flips header=0 -> 1).
    xml_query = XML_TEMPLATE.format(species=species)
    df = fetch_wbp_biomart_using_xml(xml_query)

    # 2) Normalise % identity for the scan
    df["% identity"] = pd.to_numeric(df["% identity"], errors="coerce").fillna(0.0)

    # 3) Find per-gene max and mark rows equal to that max
    mx = df.groupby("Gene stable ID")["% identity"].transform("max")
    at_max = df["% identity"].eq(mx)

    # 4) Count how many rows per gene are at the max
    tie_counts = (
        df[at_max]
        .groupby("Gene stable ID", as_index=False)
        .size()
        .sort_values("size", ascending=False)
    )

    # Keep only genes where >= 2 rows share the max (% identity ties)
    ties = tie_counts.query("size > 1")

    print(f"\n[summary] Species filter: {species}")
    print(f"[summary] Total rows fetched: {len(df):,}")
    print(f"[summary] Genes with tied max % identity: {len(ties):,}\n")

    if ties.empty:
        print("No ties found. The old reducer will appear deterministic on today's data.")
        return

    # 5) Show a sample of tie groups and their candidate rows
    print("[sample] First few genes with ties (gene_id, #rows at max):")
    print(ties.head(sample_n).to_string(index=False))

    sample_genes = set(ties.head(sample_n)["Gene stable ID"])
    sample_rows = (
        df[df["Gene stable ID"].isin(sample_genes)]
        .merge(ties, on="Gene stable ID", how="inner")
        .sort_values(["size", "Gene stable ID", "% identity"], ascending=[False, True, False])
    )

    # Keep only rows that are at the gene's max % identity
    sample_rows = sample_rows[sample_rows["% identity"].eq(
        sample_rows.groupby("Gene stable ID")["% identity"].transform("max")
    )]

    cols = ["Gene stable ID", "Human gene name", "% identity", "Homology type"]
    print("\n[sample] Candidate rows at the max for those genes:")
    print(sample_rows[cols].head(sample_n * 6).to_string(index=False))  # show a few

    # 6) Optional: write a CSV for inspection
    out = "ties_scan_sample.csv"
    sample_rows[cols + ["size"]].to_csv(out, index=False)
    print(f"\n[wrote] {out}  (inspect to see which rows compete within each gene)")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Scan live BioMart output for per-gene max %identity ties.")
    ap.add_argument("--species", default="trtricprjeb535,wubancprjna275548",
                    help="species_id_1010 value, comma-separated (old pipeline default)")
    ap.add_argument("--sample-n", type=int, default=10, help="how many tie-genes to preview")
    args = ap.parse_args()
    main(args.species, args.sample_n)
