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
is_enzyme.to_csv('pipeline/is_enzyme.csv', index=False)
overall = pd.merge(overall, is_enzyme, on='Gene stable ID', how='left')

is_kinase = testInterProGeneOntology(genomes, "GO:0004672", "is_kinase")
is_kinase.to_csv('pipeline/is_kinase.csv', index=False)
overall = pd.merge(overall, is_kinase, on='Gene stable ID', how='left')

is_ion_channel = testInterProGeneOntology(genomes, "GO:0015267", "is_ion_channel")
is_ion_channel.to_csv('pipeline/is_ion_channel.csv', index=False)
overall = pd.merge(overall, is_ion_channel, on='Gene stable ID', how='left')

is_gpcr = testInterProGeneOntology(genomes, "GO:0004930", "is_gpcr")
is_gpcr.to_csv('pipeline/is_gpcr.csv', index=False)
overall = pd.merge(overall, is_gpcr, on='Gene stable ID', how='left')

is_nuclear_receptor = testInterProGeneOntology(genomes, "GO:0004879", "is_nuclear_receptor")
is_nuclear_receptor.to_csv('pipeline/is_nuclear_receptor.csv', index=False)
overall = pd.merge(overall, is_nuclear_receptor, on='Gene stable ID', how='left')

overall.to_csv('pipeline/overall.csv', index=False)

print("Number of genes that encode enzymes:", overall[is_enzyme['is_enzyme'] == True].shape[0])
print("Number of genes that encode kinases:", overall[is_kinase['is_kinase'] == True].shape[0])
print("Number of genes that encode ion channels:", overall[is_ion_channel['is_ion_channel'] == True].shape[0])
print("Number of genes that encode gpcrs:", overall[is_gpcr['is_gpcr'] == True].shape[0])
print("Number of genes that encode nuclear receptors:", overall[is_nuclear_receptor['is_nuclear_receptor'] == True].shape[0])







