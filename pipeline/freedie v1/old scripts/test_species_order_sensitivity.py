import pandas as pd
from queryWbpBiomart import fetch_wbp_biomart_using_xml  # old client (no sort)

XML = """<?xml version="1.0" encoding="UTF-8"?>
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

def fetch(species_csv: str) -> pd.DataFrame:
    df = fetch_wbp_biomart_using_xml(XML.format(species=species_csv))
    df["% identity"] = pd.to_numeric(df["% identity"], errors="coerce").fillna(0.0)
    return df

def idxmax_reduce(dfx: pd.DataFrame) -> pd.DataFrame:
    # same as old reducer logic
    winners_idx = dfx.groupby("Gene stable ID")["% identity"].idxmax()
    return (dfx.loc[winners_idx, ["Gene stable ID", "Human gene name", "% identity"]]
              .rename(columns={"Human gene name": "winner"}))

def find_tie_genes(dfx: pd.DataFrame) -> set:
    mx = dfx.groupby("Gene stable ID")["% identity"].transform("max")
    at_max = dfx["% identity"].eq(mx)
    tie_counts = (dfx[at_max].groupby("Gene stable ID").size())
    return set(tie_counts[tie_counts > 1].index)

def main():
    A_then_B = "trtricprjeb535,wubancprjna275548"
    B_then_A = "wubancprjna275548,trtricprjeb535"

    df_A = fetch(A_then_B)
    df_B = fetch(B_then_A)

    # Reduce winners with the old idxmax approach
    wA = idxmax_reduce(df_A)
    wB = idxmax_reduce(df_B)

    # Compare winners where % identity is identical (focus on true ties)
    merged = wA.merge(wB, on=["Gene stable ID", "% identity"], suffixes=("_A","_B"))
    changed = merged[merged["winner_A"] != merged["winner_B"]]

    # Optional: restrict to genes that actually have ties in the data
    ties = find_tie_genes(df_A)
    changed_ties = changed[changed["Gene stable ID"].isin(ties)]

    print(f"[summary] A,B vs B,A total genes compared: {len(merged):,}")
    print(f"[summary] Genes with different winners (same % identity): {len(changed):,}")
    print(f"[summary] …of which are known tie-genes: {len(changed_ties):,}\n")

    if not changed_ties.empty:
        print("[sample] Changed winners (first 20):")
        print(changed_ties.head(20).to_string(index=False))

        # Evidence rows for the first few changed genes (what were the candidates?)
        sample_ids = set(changed_ties["Gene stable ID"].head(5))
        ev = (df_A[df_A["Gene stable ID"].isin(sample_ids)]
              .copy()
              .sort_values(["Gene stable ID", "% identity"], ascending=[True, False]))
        ev_cols = ["Gene stable ID", "Human gene name", "% identity", "Homology type"]
        ev.to_csv("species_order_changes_evidence.csv", index=False, columns=ev_cols)
        changed_ties.to_csv("species_order_changed_winners.csv", index=False)
        print("\n[wrote] species_order_changed_winners.csv")
        print("[wrote] species_order_changes_evidence.csv")
    else:
        print("No changed winners detected. There are ties, but species order didn't alter row order today.")

if __name__ == "__main__":
    main()
