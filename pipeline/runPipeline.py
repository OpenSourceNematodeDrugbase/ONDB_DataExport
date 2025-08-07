import pandas as pd
from populateGeneList import *

geneList = retrieveGeneListFromWbpBiomart()

print(geneList)

geneList.to_csv('pipeline/geneList.csv', index=False)
