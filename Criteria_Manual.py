import pandas as pd
import uuid
import re
import json
import os

def sanitize_key(key):
    key = key.lower()
    key = re.sub(r'[.#$/\[\]\s%()]+', '_', key)
    return key.strip('_')

# ✅ Correct base path
base_path = r"F:\ONDB_DATA"

# ✅ File names directly in base_path (no 'Data' subfolder)
file_names = [
    "mart_export.tsv",
    "mart_export_paralogue.tsv",
    "mart_export_GoTerm.tsv"
]

# ✅ Build full paths correctly
file_paths = [os.path.join(base_path, name) for name in file_names]

# ✅ Optional: check for file existence
for path in file_paths:
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

# ✅ Load and combine data
df_all = pd.concat([pd.read_csv(fp, sep='\t') for fp in file_paths], ignore_index=True)
print(f"Loaded {len(df_all)} total records from {len(file_paths)} files")

# ✅ Assign unique IDs
df_all['id'] = [str(uuid.uuid4()) for _ in range(len(df_all))]

# ✅ Convert to dicts
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


def check_data_match(column_names, criteria_key, check_input):
    """
    Check if any of the values in check_input are present in the given columns.
    Sets the criteria_key to True if a match is found, else False.

    Args:
        column_names (list of str): Columns to check for matches.
        criteria_key (str): The key to assign in records indicating match status.
        check_input (list or set): Values to look for in the columns.
    """
    existing_columns = [col for col in column_names if column_exists(col)]
    if not existing_columns:
        print(f"None of the columns for '{criteria_key}' exist. Skipping...")
        return

    # Normalize check_input into a set for efficient lookup
    check_set = set(check_input)

    # Track ids that have matches
    matched_ids = set()

    for record in records:
        # For each record, check if any of the specified columns contain any of the check_input values
        for col in existing_columns:
            value = record.get(col)
            if pd.isnull(value):
                continue

            # Sometimes GO terms may be concatenated with delimiters; split if necessary
            terms = re.split(r'[;, ]+', str(value))
            if any(term in check_set for term in terms):
                matched_ids.add(record["id"])
                break  # No need to check further columns for this record

    # Assign True if matched, False otherwise
    for record in records:
        record[criteria_key] = (record["id"] in matched_ids)

    # ✅ Print how many matched
    print(f"Number of records with {criteria_key}: {len(matched_ids)}")


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

go_terms_to_check = ["GO:0002168", "GO:0002119", "GO:0061062", "GO:0007275"]
check_data_match(["GO term accession"], "linked_to_larval_development", go_terms_to_check)

# --- Export ---
export_data_json("records.json")