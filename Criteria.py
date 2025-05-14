import uuid
import pandas as pd
import requests
from io import StringIO
import re

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

df_response = pd.read_csv(StringIO(response.text), sep="\t")

# Assuming df_response is your DataFrame
print(f"Querying {len(df_response)} records...")

# Generate a unique 'id' for each row in the DataFrame
df_response['id'] = [str(uuid.uuid4()) for _ in range(len(df_response))]

# Convert the DataFrame to a list of dictionaries (records)
records = df_response.to_dict(orient="records")

def look_up_keyword(df, column_name, keyword):
    results = df[df[column_name].str.contains(keyword, case=False, na=False)]
    return results

def look_up_null(df, column_name):
    results = []

    # Iterate over rows where the specified column has NaN values
    for row_index, row in df.iterrows():
        if pd.isnull(row[column_name]):
            results.append((row_index, row))  # Store index and entire row as a tuple

    return results

def look_up_value(df, greater_than, threshold, column_name):
    results = []

    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in DataFrame.")

    for row_index, row in df.iterrows():
        value = row[column_name]
        if pd.notnull(value):  # Skip Null Values
            if greater_than and value > threshold:
                results.append((row_index, column_name))
            elif not greater_than and value < threshold:
                results.append((row_index, column_name))
    return results

#This is an example of a criteria function
def human_simularity():
    null_human_orthologue = look_up_null(df_response, "Human protein stable ID")

    # Extract the 'id' from the null_human_orthologue list of tuples
    null_ids = {row[1]["id"] for row in null_human_orthologue}

    print(f"Number of drug targets that do not have a similar protein to humans: {len(null_human_orthologue)}")

    # Iterate through the records and update 'is_in_tuple' based on whether the ID is in null_ids
    for record in records:
        if record["id"] in null_ids:
            record["orthologue_comparison"] = False
        else:
            record["orthologue_comparison"] = True



human_simularity()

















