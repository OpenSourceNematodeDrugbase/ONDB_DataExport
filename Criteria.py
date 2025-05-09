import pandas as pd
import os
import requests
from io import StringIO
import re

def look_up_keyword(df, column_name, keyword):
    results = df[df[column_name].str.contains(keyword, case=False, na=False)]
    return results

def is_data_null(df, column_name):
    for index, row in df.iterrows():
        if row[column_name].isnull():
            return True
        else:
            return False



# Wormbase Parasite BioMart URL
biomart_url = "https://parasite.wormbase.org/biomart/martservice?query="

xml_query = """
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "wbps_gene" interface = "default" >
		<Filter name = "species_id_1010" value = "wubancprjna275548"/>
		<Attribute name = "production_name_1010" />
		<Attribute name = "hsapiens_gene" />
		<Attribute name = "hsapiens_gene_name" />
		<Attribute name = "hsapiens_homolog_ensembl_peptide" />
		<Attribute name = "hsapiens_orthology_type" />
		<Attribute name = "hsapiens_homolog_perc_id" />
		<Attribute name = "hsapiens_homolog_perc_id_r1" />
		<Attribute name = "gene_biotype" />
		<Attribute name = "wbps_gene_id" />
		<Attribute name = "caelegprjna13758_gene" />
		<Attribute name = "caelegprjna13758_gene_name" />
		<Attribute name = "caelegprjna13758_homolog_perc_id" />
	</Dataset>
</Query>
"""

# remove newlines in the xml_query string
xml_query = xml_query.replace('\n', '').replace('\r', '')

# change header = "0" to header = "1" to get the column names in the first row of the biomart output

xml_query = re.sub(r'header\s*=\s*"0"', 'header="1"', xml_query)
response = requests.get(biomart_url + xml_query)

df = pd.read_csv(StringIO(response.text), sep="\t")

length = 0
nullLength = 0
validLength = 0

for index, row in df.iterrows():
    is_nan = pd.isnull(df.loc[index, "Human gene stable ID"])
    print(f"Gene: {df.loc[index, "Gene stable ID"]}, Column: Human gene stable ID, Is Null: {is_nan}")

    length += 1

    if is_nan:
        nullLength += 1
    else:
        validLength += 1

print(f"Total Length: {length}")
print(f"Total Null Length: {nullLength}")
print(f"Total Valid Length: {validLength}")

gene = "wuchereria_bancrofti_prjna275548"
prop = "Human gene stable ID"
val = df.loc[0, prop]

proteins = look_up_keyword(df, "Gene biotype", "protein_coding")

print(f"Results for Proteins: {len(proteins)}")

identityValid = 0

for index, row in df.iterrows():
    if row["% identity"] < 10:
        print(f"Gene: {df.loc[index, "Gene stable ID"]}, Column: %Identity: {row["% identity"]}")
        identityValid += 1

print(f"Total Valid Identity: {identityValid}")

#Determine Score
#Save Attributes to Gene Output
#For each class, export the data and scores

#Class (Gene)
#identity%
#Gene biotype













