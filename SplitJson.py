import json
import os
import argparse

def split_json(input_file, output_dir, chunk_size):
    # Load the entire JSON file
    with open(input_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    if not isinstance(data, list):
        raise ValueError("Expected a JSON array at the root level.")

    os.makedirs(output_dir, exist_ok=True)
    total_chunks = (len(data) + chunk_size - 1) // chunk_size

    for i in range(total_chunks):
        chunk = data[i * chunk_size : (i + 1) * chunk_size]
        output_path = os.path.join(output_dir, f"output_{i + 1}.json")
        with open(output_path, 'w', encoding='utf-8') as out_f:
            json.dump(chunk, out_f, ensure_ascii=False, indent=2)
        print(f"Wrote {len(chunk)} records to {output_path}")

    print(f"\n✅ Successfully split into {total_chunks} file(s) in '{output_dir}'.")

split_json("records.json", "JsonExport", 25000)
