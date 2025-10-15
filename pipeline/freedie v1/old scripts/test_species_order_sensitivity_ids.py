import pandas as pd
from queryWbpBiomart import fetch_wbp_biomart_using_xml  # old client

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
    df = fetch_wbp_biomart_using_xml(XML.format(species=species_csv)).reset_index(drop=True)
    df["row_index"] = df.index  # preserve true input order
    df["% identity"] = pd.to_numeric(df["% identity"], errors="coerce").fillna(0.0)
    return df

def first_at_max_per_gene(dfx: pd.DataFrame) -> pd.DataFrame:
    # rows at max % identity within each gene
    mx = dfx.groupby("Gene stable ID")["% identity"].transform("max")
    at_max = dfx[dfx["% identity"].eq(mx)]
    # pick the first-at-max by the original row order
    first = (at_max.sort_values(["Gene stable ID", "row_index"])
                  .drop_duplicates(subset=["Gene stable ID"], keep="first"))
    # robust winner label (prefer stable ID)
    first["winner_id"] = first["Human gene stable ID"].fillna("").astype(str)
    first["winner_name"] = first["Human gene name"].fillna("").astype(str)
    first["winner"] = first["winner_id"].where(first["winner_id"] != "", first["winner_name"])
    return first[["Gene stable ID", "% identity", "winner", "winner_id", "winner_name", "row_index"]]

def main():
    A_then_B = "trtricprjeb535,wubancprjna275548,trregeprjeb44434"

    B_then_A = "trregeprjeb44434,wubancprjna275548,trtricprjeb535"


    dfA = fetch(A_then_B)
    dfB = fetch(B_then_A)

    wA = first_at_max_per_gene(dfA)
    wB = first_at_max_per_gene(dfB)

    merged = wA.merge(wB, on=["Gene stable ID", "% identity"], suffixes=("_A","_B"))
    changed = merged[merged["winner_A"].fillna("") != merged["winner_B"].fillna("")]

    print(f"[summary] compared genes: {len(merged):,}")
    print(f"[summary] different winners (same % identity): {len(changed):,}")

    if not changed.empty:
        print("[sample] changed winners (first 20):")
        print(changed[["Gene stable ID","% identity","winner_A","winner_B","row_index_A","row_index_B"]]
              .head(20).to_string(index=False))

        # Evidence: candidate rows at max, with their row_index
        ids = set(changed["Gene stable ID"].head(10))
        evA = dfA[dfA["Gene stable ID"].isin(ids)]
        evB = dfB[dfB["Gene stable ID"].isin(ids)]
        # keep only candidates at max
        evA = evA[evA["% identity"].eq(evA.groupby("Gene stable ID")["% identity"].transform("max"))]
        evB = evB[evB["% identity"].eq(evB.groupby("Gene stable ID")["% identity"].transform("max"))]
        cols = ["Gene stable ID","Human gene stable ID","Human gene name","% identity","Homology type","row_index"]
        evA[cols].to_csv("species_order_changes_evidence_A.csv", index=False)
        evB[cols].to_csv("species_order_changes_evidence_B.csv", index=False)
        changed.to_csv("species_order_changed_winners.csv", index=False)
        print("\n[wrote] species_order_changed_winners.csv")
        print("[wrote] species_order_changes_evidence_A.csv")
        print("[wrote] species_order_changes_evidence_B.csv")
    else:
        print("No changed winners: server returned identical per-gene row order for A,B and B,A today (ties still exist).")

if __name__ == "__main__":
    main()


# [summary] compared genes: 23,804
# [summary] different winners (same % identity): 0
# No changed winners: server returned identical per-gene row order for A,B and B,A today (ties still exist).