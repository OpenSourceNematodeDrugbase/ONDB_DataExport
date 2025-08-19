import pandas as pd
import polars as pl
import numpy as np

from populateGeneList import *
from wbpHumanOrthologues import *
from testInterProGeneOntology import *
from testInterProStringSearch import *

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

# search for GPCR or G-protein coupled receptor, but not GPCR kinase  in the InterPro domain annotations
# regular expressions for polars can't have lookaround hence basic approach below that needs to be checked carefully
is_gpcr = testInterProStringSearch(genomes, 'G-protein coupled receptor|GPCR.[^k]', "is_gpcr")
is_gpcr = is_gpcr.to_pandas()
is_gpcr.to_csv('pipeline/is_gpcr.csv', index=False)
overall = pd.merge(overall, is_gpcr, on='Gene stable ID', how='left')

is_nuclear_receptor = testInterProGeneOntology(genomes, "GO:0004879", "is_nuclear_receptor")
is_nuclear_receptor.to_csv('pipeline/is_nuclear_receptor.csv', index=False)
overall = pd.merge(overall, is_nuclear_receptor, on='Gene stable ID', how='left')

# is privileged target family if any of these are true
overall['is_privileged_target_family'] = overall[['is_enzyme', 'is_kinase', 'is_ion_channel', 'is_gpcr', 'is_nuclear_receptor']].any(axis=1)
overall['is_privileged_target_family_evidence'] = np.where(overall.is_privileged_target_family==True, 'Meets the following criteria: ' + overall[['is_enzyme', 'is_kinase', 'is_ion_channel', 'is_gpcr', 'is_nuclear_receptor']].apply(
    lambda row: ", ".join(row.index[row]), axis=1
), "Does not meet the following criteria:  'is_enzyme', 'is_kinase', 'is_ion_channel', 'is_gpcr', 'is_nuclear_receptor'")


overall.to_csv('pipeline/overall.csv', index=False)

print("Number of genes that encode enzymes:", overall[is_enzyme['is_enzyme'] == True].shape[0])
print("Number of genes that encode kinases:", overall[is_kinase['is_kinase'] == True].shape[0])
print("Number of genes that encode ion channels:", overall[is_ion_channel['is_ion_channel'] == True].shape[0])
print("Number of genes that encode gpcrs:", overall[is_gpcr['is_gpcr'] == True].shape[0])
print("Number of genes that encode nuclear receptors:", overall[is_nuclear_receptor['is_nuclear_receptor'] == True].shape[0])

