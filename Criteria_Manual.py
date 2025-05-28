import pandas as pd
import uuid
import re
import json
import os

def sanitize_key(key):
    key = key.lower()
    key = re.sub(r'[.#$/\[\]\s%()]+', '_', key)
    return key.strip('_')

# --- Load and Combine Multiple Files ---
file_paths = [
    "Data/mart_export.tsv",
    "Data/mart_export_paralogue.tsv"
]  # Add all files you want to combine here

df_all = pd.concat([pd.read_csv(fp, sep='\t') for fp in file_paths], ignore_index=True)
print(f"Loaded {len(df_all)} total records from {len(file_paths)} files")

# --- Assign ID Once Per Unique Record ---
df_all['id'] = [str(uuid.uuid4()) for _ in range(len(df_all))]

# --- Convert to Dicts ---
records = df_all.to_dict(orient="records")

def column_exists(column_name):
    return column_name in df_all.columns

def look_up_null(df, column_name, return_ids_only=False):
    results = []
    for row_index, row in df.iterrows():
        if pd.isnull(row[column_name]):
            if return_ids_only:
                results.append(row["id"])
            else:
                results.append((row_index, row))
    return results

def check_data_assigned(column_names, criteria_key, null_is_true):
    existing_columns = [col for col in column_names if column_exists(col)]
    if not existing_columns:
        print(f"None of the columns for '{criteria_key}' exist. Skipping...")
        return

    null_ids = set()
    for column_name in existing_columns:
        null_ids.update(look_up_null(df_all, column_name, True))

    print(f"Number of drug targets that DO NOT have {criteria_key} assigned: {len(null_ids)}")

    for record in records:
        record[criteria_key] = (record["id"] in null_ids) == null_is_true

def export_data_json(filename):
    clean_records = []

    for record in records:
        clean_record = {"id": sanitize_key(record["id"])}
        for key, value in record.items():
            if key == "id":
                continue
            if pd.isnull(value):
                value = "NaN"
            clean_record[sanitize_key(key)] = value
        clean_records.append(clean_record)

    with open(filename, 'w', encoding='utf-8') as file:
        json.dump(clean_records, file, indent=4)
    print(f"Data exported to {filename}")

# --- Run Checks ---
check_data_assigned(["Human protein stable ID"], "similar_protein_in_humans", null_is_true=False)
check_data_assigned(["InterPro ID"], "has_known_protein_domain", null_is_true=False)
check_data_assigned(["GO term accession", "GO term name", "GO term evidence code"], "has_gene_ontology_functional_annotation", null_is_true=False)
check_data_assigned(["Paralogue gene stable ID"], "has_paralogue_id", null_is_true=False)

# --- Export ---
export_data_json("records.json")