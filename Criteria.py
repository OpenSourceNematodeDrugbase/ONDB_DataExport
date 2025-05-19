import uuid
import pandas as pd
import requests
from io import StringIO
import re
import json
import sys

# Wormbase Parasite BioMart URL
biomart_url = "https://parasite.wormbase.org/biomart/martservice?query="

xml_query = """
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "wbps_gene" interface = "default" >
        <Filter name = "species_id_1010" value = "wubancprjna275548"/>
        <Attribute name = "production_name_1010" />
        <Attribute name = "wbps_gene_id" />
        <Attribute name = "wbps_paralog_gene" />
        <Attribute name = "wbps_gene__paralog__dm_perc_id_4015" />
        <Attribute name = "wbps_gene__paralog__dm_perc_id_4015_r1" />
    </Dataset>
</Query>
"""

# remove newlines in the xml_query string
xml_query = xml_query.replace('\n', '').replace('\r', '')

# change header = "0" to header = "1" to get the column names in the first row of the biomart output

xml_query = re.sub(r'header\s*=\s*"0"', 'header="1"', xml_query)
response = requests.get(biomart_url + xml_query, stream=True)

total_size = int(response.headers.get('content-length', 0))
chunk_size = 1024
downloaded = 0
data = ""

print("Downloading BioMart data...")

for chunk in response.iter_content(chunk_size=chunk_size):
        if chunk:
            downloaded += len(chunk)
            data += chunk.decode("utf-8")
            sys.stdout.write(f"\rDownloaded {downloaded / 1024:.1f} KB")
            sys.stdout.flush()

print("\nDownload complete.")

df_response = pd.read_csv(StringIO(data), sep="\t")

# Assuming df_response is your DataFrame
print(f"Querying {len(df_response)} records...")

# Generate a unique 'id' for each row in the DataFrame
df_response['id'] = [str(uuid.uuid4()) for _ in range(len(df_response))]

# Convert the DataFrame to a list of dictionaries (records)
records = df_response.to_dict(orient="records")

def column_exists(column_name):
    return column_name in df_response.columns

def look_up_keyword(df, column_name, keyword):
    results = df[df[column_name].str.contains(keyword, case=False, na=False)]
    return results

def look_up_null(df, column_name, return_ids_only=False):
    results = []

    for row_index, row in df.iterrows():
        if pd.isnull(row[column_name]):
            if return_ids_only:
                results.append(row["id"])
            else:
                results.append((row_index, row))

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

def check_data_assigned(column_names, criteria_key, null_is_true):
    null_ids = set()

    for column_name in column_names:
        if column_exists(column_name):
            null_ids.update(look_up_null(df_response, column_name, True))
        else:
            print(f"Column '{column_name}' does not exist.")
            null_is_true = True

    print(f"Number of drug targets that DO NOT have {criteria_key} assigned: {len(null_ids)}")

    for record in records:
        if record["id"] in null_ids:
            record[criteria_key] = null_is_true
        else:
            record[criteria_key] = not null_is_true


def sanitize_key(key):
    # Lowercase for consistency
    key = key.lower()

    # Replace invalid Firebase chars (including space, %, (), etc.) with underscore
    key = re.sub(r'[.#$/\[\]\s%()]+', '_', key)

    # Remove leading/trailing underscores
    key = key.strip('_')

    return key

def export_data_json(filename):
    clean_records = []

    for record in records:
        clean_record = {}
        for key, value in record.items():
            # Skip 'id' field
            if key == 'id':
                continue

            # Replace NaN with "NaN" string
            if pd.isnull(value):
                value = "NaN"

            # Sanitize the key for Firebase
            safe_key = sanitize_key(key)

            clean_record[safe_key] = value

        clean_records.append(clean_record)

    # Write to JSON
    with open(filename, 'w', encoding='utf-8') as file:
        json.dump(clean_records, file, indent=4)

    print(f"Data exported to {filename}")

print(df_response.columns)

check_data_assigned(["Human protein stable ID"], "similar_protein_in_humans", null_is_true=False) #Human Similarity
check_data_assigned(["InterPro ID"], "has_known_protein_domain", null_is_true=False) #Has InterPro ID
check_data_assigned(["GO term accession", "GO term name", "GO term evidence code"], "has_gene_ontology_functional_annotation", null_is_true=False) #Gene Ontology (GO) functional annotation
check_data_assigned(["Paralogue gene stable ID"], "has_gene_ontology_functional_annotation", null_is_true=False) #Has paralogue ID

export_data_json("records.json")



















