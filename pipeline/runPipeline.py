import pandas as pd
from populateGeneList import *
from wbpHumanOrthologues import *

# define list of genomes to query

# for now, we will just use the most recent Wuchereria bancrofti genome
genomes = "wubancprjna275548" 


geneList = retrieveGeneListFromWbpBiomart(genomes)
geneList.to_csv('pipeline/geneList.csv', index=False)
print("Number of genes retrieved:", geneList.shape[0])


wbpHumanOrthologues = queryWbpHumanOrthologues(genomes)
wbpHumanOrthologues.to_csv('pipeline/wbpHumanOrthologues.csv', index=False)
print("Number of genes lacking WBP human orthologue:", wbpHumanOrthologues[wbpHumanOrthologues['lacks_WBP_human_orthologue'] == True].shape[0])
print("Number of genes lacking best WBP human orthologue >= 40% identity:", wbpHumanOrthologues[wbpHumanOrthologues['best_WBP_human_orthologue_lt_40pct_identity'] == True].shape[0])

overall = pd.merge(geneList, wbpHumanOrthologues, on='Gene stable ID', how='left')
overall.to_csv('pipeline/overall.csv', index=False)
