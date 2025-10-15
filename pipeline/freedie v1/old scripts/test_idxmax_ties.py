import pandas as pd

def old_reduce_human(df):
    df = df.copy()
    df["% identity"] = df["% identity"].fillna(0)
    return df.loc[df.groupby("Gene stable ID")["% identity"].idxmax()][
        ["Gene stable ID","Human gene name","% identity"]
    ]

df1 = pd.DataFrame({
    "Gene stable ID": ["WBGene0001","WBGene0001"],
    "Human gene name": ["HUMAN_A","HUMAN_B"],
    "% identity": [55.0,55.0],  # tie
})
print("Order #1:\n", old_reduce_human(df1), "\n")

df2 = df1.iloc[::-1].reset_index(drop=True)  # reversed input
print("Order #2:\n", old_reduce_human(df2))
