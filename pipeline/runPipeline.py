import pandas as pd
from populateGeneList import *
from wbpHumanOrthologues import *
from testInterProGeneOntology import *

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

is_enzyme = testInterProGeneOntology(genomes, "GO:0003824", "is_enzyme")
is_enzyme.to_csv('pipeline/is_Enzyme.csv', index=False)
print("Number of genes that are enzymes:", is_enzyme[is_enzyme['is_enzyme'] == True].shape[0])


overall = pd.merge(overall, is_enzyme, on='Gene stable ID', how='left')
overall.to_csv('pipeline/overall.csv', index=False)




