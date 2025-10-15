import pandas as pd

# Two human orthologues tie on % identity for the same parasite gene
df = pd.DataFrame({
    "Gene stable ID": ["G1", "G1"],
    "Human gene name": ["HUMAN_B", "HUMAN_A"],
    "% identity": [55.0, 55.0],  # tie
})

def pick_idxmax(dfx: pd.DataFrame) -> pd.DataFrame:
    dfx = dfx.copy()
    dfx["% identity"] = dfx["% identity"].fillna(0)
    # This mirrors the old reducer logic in wbpHumanOrthologues (legacy):
    # df.loc[df.groupby('Gene stable ID')['% identity'].idxmax()]
    return dfx.loc[dfx.groupby("Gene stable ID")["% identity"].idxmax()][
        ["Gene stable ID", "Human gene name", "% identity"]
    ]

print("Original order winner:\n", pick_idxmax(df), "\n")
print("Reversed order winner:\n", pick_idxmax(df.iloc[::-1].reset_index(drop=True)))

#When two (or more) rows tie on % identity within a gene, the choice depends on whatever row arrives first.