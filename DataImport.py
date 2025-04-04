import pandas as pd
import os

def look_up_keyword(df, column_name, keyword):
    results = df[df[column_name].str.contains(keyword, case=False, na=False)]
    print(f"Results for lethal: {len(results)}")
    return results

def export_json(results, output_file_name):
    # Gets the root directory of the script
    root_directory = os.getcwd()
    output_path = os.path.join(root_directory, output_file_name)

    # Export results to JSON in root directory
    results.to_json(output_path, orient="records", indent=4)

    print(f"Results exported to {output_path}")


mart_df = pd.read_csv('mart_export.txt')
print(mart_df.head())

sm_df = pd.read_csv("simplemine_results.txt", sep="\t")
print(sm_df.head())

lethalResults = look_up_keyword(sm_df, "RNAi Phenotype Observed", 'lethal')
print(lethalResults)

export_json(lethalResults, "lethalResults.json")




