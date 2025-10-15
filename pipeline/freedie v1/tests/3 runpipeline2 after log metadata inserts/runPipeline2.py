#coded on the 16 sept - links the alphafold script - changed to a defined species once using upper case GENOMES
# works but needs to sort the JSON metadata within logs folder that is not being created

import os, sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)


import pandas as pd
import polars as pl
import numpy as np

from populateGeneList import *
from wbpHumanOrthologues import *
from testInterProGeneOntology import *
from testInterProStringSearch import *
from alphafoldAvailability import queryAlphaFoldAvailability

# define list of genomes to query

# for now, we will just use the most recent Wuchereria bancrofti genome and the Trichuris trichiura genome
GENOMES = "trtricprjeb535,wubancprjna275548" 

geneList = retrieveGeneListFromWbpBiomart(GENOMES)


geneList.to_csv('geneList.csv', index=False)
print("Number of genes retrieved:", geneList.shape[0])




wbpHumanOrthologues = queryWbpHumanOrthologues(GENOMES)
wbpHumanOrthologues.to_csv('wbpHumanOrthologues.csv', index=False)
print("Number of genes lacking WBP human orthologue:", wbpHumanOrthologues[wbpHumanOrthologues['lacks_WBP_human_orthologue'] == True].shape[0])
print("Number of genes lacking best WBP human orthologue >= 40% identity:", wbpHumanOrthologues[wbpHumanOrthologues['best_WBP_human_orthologue_lt_40pct_identity'] == True].shape[0])

overall = pd.merge(geneList, wbpHumanOrthologues, on='Gene stable ID', how='left')



is_enzyme = testInterProGeneOntology(GENOMES, "GO:0003824", "is_enzyme")
is_enzyme.to_csv('is_enzyme.csv', index=False)
overall = pd.merge(overall, is_enzyme, on='Gene stable ID', how='left')

is_kinase = testInterProGeneOntology(GENOMES, "GO:0004672", "is_kinase")
is_kinase.to_csv('is_kinase.csv', index=False)
overall = pd.merge(overall, is_kinase, on='Gene stable ID', how='left')

is_ion_channel = testInterProGeneOntology(GENOMES, "GO:0015267", "is_ion_channel")
is_ion_channel.to_csv('is_ion_channel.csv', index=False)
overall = pd.merge(overall, is_ion_channel, on='Gene stable ID', how='left')

# search for GPCR or G-protein coupled receptor, but not GPCR kinase  in the InterPro domain annotations
# regular expressions for polars can't have lookaround hence basic approach below that needs to be checked carefully
is_gpcr = testInterProStringSearch(GENOMES, 'G-protein coupled receptor|GPCR.[^k]', "is_gpcr")
is_gpcr = is_gpcr.to_pandas()
is_gpcr.to_csv('is_gpcr.csv', index=False)
overall = pd.merge(overall, is_gpcr, on='Gene stable ID', how='left')

is_nuclear_receptor = testInterProGeneOntology(GENOMES, "GO:0004879", "is_nuclear_receptor")
is_nuclear_receptor.to_csv('is_nuclear_receptor.csv', index=False)
overall = pd.merge(overall, is_nuclear_receptor, on='Gene stable ID', how='left')

# is privileged target family if any of these are true
overall['is_privileged_target_family'] = overall[['is_enzyme', 'is_kinase', 'is_ion_channel', 'is_gpcr', 'is_nuclear_receptor']].any(axis=1)
overall['is_privileged_target_family_evidence'] = np.where(overall.is_privileged_target_family==True, 'Meets the following criteria: ' + overall[['is_enzyme', 'is_kinase', 'is_ion_channel', 'is_gpcr', 'is_nuclear_receptor']].apply(
    lambda row: ", ".join(row.index[row]), axis=1
), "Does not meet the following criteria:  'is_enzyme', 'is_kinase', 'is_ion_channel', 'is_gpcr', 'is_nuclear_receptor'")

overall.to_csv('overall.csv', index=False)


import importlib, sys, os
print("[debug] CWD=", os.getcwd())
spec = importlib.util.find_spec("utils_logging")
print("[debug] utils_logging spec:", spec)
if spec is not None:
    mod = importlib.import_module("utils_logging")
    print("[debug] utils_logging loaded from:", mod.__file__)



#### NEW INSERT LOGGIN METADATA ####
# ---- AlphaFold availability (logs into pipeline/logs/RUN_METADATA.json) ----
alphafold_df = queryAlphaFoldAvailability(GENOMES, log=True)  # uses the same species list

# Keep a CSV for auditing (the logger already records this path too)
alphafold_df.to_csv('alphafold_availability.csv', index=False)

# Merge onto your main table
overall = pd.merge(
    overall,
    alphafold_df[['Gene stable ID', 'alphafold_available']],
    on='Gene stable ID',
    how='left'
).copy()

overall['alphafold_available'] = overall['alphafold_available'].fillna(False)

print("Number of genes WITH AlphaFold:", int((overall['alphafold_available'] == True).sum()))
print("Number of genes WITHOUT AlphaFold:", int((overall['alphafold_available'] == False).sum()))


export = overall[['Gene stable ID', 'Genome name', 'Transcript stable ID', 'Gene description',
                   'alphafold_available',                                 # <-- NEW
                   'lacks_WBP_human_orthologue', 'lacks_WBP_human_orthologue_evidence',
                   'best_WBP_human_orthologue_lt_40pct_identity', 'best_WBP_human_orthologue_lt_40pct_identity_evidence',
                   'is_enzyme', 'is_enzyme_evidence', 
                   'is_kinase', 'is_kinase_evidence', 
                   'is_ion_channel', 'is_ion_channel_evidence', 
                   'is_gpcr', 'is_gpcr_evidence', 
                   'is_nuclear_receptor', 'is_nuclear_receptor_evidence',
                   'is_privileged_target_family', 'is_privileged_target_family_evidence']]


export.to_csv('export.csv', index=False)
export.to_json('export.json', orient='records', indent=4)

print("Number of genes that encode enzymes:", overall[is_enzyme['is_enzyme'] == True].shape[0])
print("Number of genes that encode kinases:", overall[is_kinase['is_kinase'] == True].shape[0])
print("Number of genes that encode ion channels:", overall[is_ion_channel['is_ion_channel'] == True].shape[0])
print("Number of genes that encode gpcrs:", overall[is_gpcr['is_gpcr'] == True].shape[0])
print("Number of genes that encode nuclear receptors:", overall[is_nuclear_receptor['is_nuclear_receptor'] == True].shape[0])

