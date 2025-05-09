import pandas as pd
import os

def look_up_keyword(df, column_name, keyword):
    results = df[df[column_name].str.contains(keyword, case=False, na=False)]
    print(f"Results for {keyword}: {len(results)}")
    return results

def export_json(results, output_file_name):
    # Gets the root directory of the script
    root_directory = os.getcwd()
    output_path = os.path.join(root_directory, output_file_name)

    # Export results to JSON in root directory
    results.to_json(output_path, orient="records", indent=4)

    print(f"Results exported to {output_path}")

def check_if_protein_or_enzyme(row):
    is_protein = False
    is_enzyme = False

    if "protein_coding" in str(row['Biotype']).lower():
        is_protein = True

    # Check 'WormPep' - If this is not empty, it's likely a protein
    if pd.notna(row['WormPep']) and row['WormPep'].strip() != "":
        is_protein = True

    # Check 'Protein Domain' for enzymatic domains (e.g., kinase, protease)
    if pd.notna(row['Protein Domain']) and any(domain in str(row['Protein Domain']).lower() for domain in
                                               ['kinase', 'protease', 'oxidase', 'phosphatase']):
        is_enzyme = True

    # Check 'UniProt' - If there's an entry here, it's likely a protein
    if pd.notna(row['UniProt']) and row['UniProt'].strip() != "":
        is_protein = True

    # Check 'Gene Ontology Association' for enzyme-related GO terms
    if pd.notna(row['Gene Ontology Association']) and any(
            term in str(row['Gene Ontology Association']).lower() for term in
            ['catalytic', 'enzyme', 'oxidoreductase', 'hydrolase']):
        is_enzyme = True

    return is_protein, is_enzyme

mart_df = pd.read_csv('mart_export.txt')

sm_df = pd.read_csv("simplemine_results.txt", sep="\t")

#Add column 'isProtein' and 'isEnzyme' and apply values for each row
sm_df['Is_Protein'], sm_df['Is_Enzyme'] = zip(*sm_df.apply(check_if_protein_or_enzyme, axis=1))
#print(sm_df[['Is_Protein', 'Is_Enzyme']])

lethalResults = look_up_keyword(sm_df, "RNAi Phenotype Observed", 'lethal')
aceResults = look_up_keyword(sm_df, "RNAi Phenotype Observed", 'ase')

export_json(lethalResults, "lethalResults.json")




