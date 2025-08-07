import pandas as pd
from queryWbpBiomart import fetch_wbp_biomart_using_xml


def testIsEnzyme(genomes):
# first we fetch all the interpro annotations for each gene

    xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
              
        <Dataset name = "wbps_gene" interface = "default" >
        <Filter name = "species_id_1010" value = "{genomes}"/>
		<Attribute name = "wbps_gene_id" />
		<Attribute name = "interpro_id" />
		<Attribute name = "interpro_description" />
        </Dataset>
    </Query>
    """

    # Fetch the data using the XML query
    df = fetch_wbp_biomart_using_xml(xml_query)

    # now we need to find the GO terms that are related to enzymes
    # we will use the GO term "catalytic activity" (GO:0003824) as a proxy for enzymes
    # and we will use the GO terms that descend from this term to find all enzymes
    # we will use the goatools library to fetch the GO terms
    
    # fetch gene ontology go-basic.obo file if not already present
    import os
    if not os.path.exists("pipeline/go-basic.obo"):
        url = 'https://purl.obolibrary.org/obo/go/go-basic.obo'
        try:
            import urllib.request
            urllib.request.urlretrieve(url, "pipeline/go-basic.obo")
        except Exception as e:
            print(f"Error fetching GO ontology: {e}")

    # Load the ontology
    from goatools import obo_parser
    go_dag = obo_parser.GODag("pipeline/go-basic.obo")

    # Find all GO terms descending from enzymes
    root_go = "GO:0003824"
    descendants = go_dag[root_go].get_all_children()
    descendants.add(root_go)  # Include the parent itself

    # fetch interpro2go mapping file if not already present
    if not os.path.exists("pipeline/interpro2go"):
        url = 'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go'
        try:
            import urllib.request
            urllib.request.urlretrieve(url, "pipeline/interpro2go")
        except Exception as e:
            print(f"Error fetching InterPro to GO mapping: {e}")

    import re
    # Read the interpro2go mapping file
    # and create a dictionary mapping InterPro IDs to GO terms

    interpro2go = {}
    with open("pipeline/interpro2go", 'r') as file:
        for line in file:
            if line.startswith("InterPro:"):
                # Example line:
                # InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:DNA binding ; GO:0003677
                parts = line.strip().split('>')
                if len(parts) != 2:
                    continue  # skip malformed lines

                ipr_part = parts[0].strip()      # "InterPro:IPR000003 Retinoid X receptor/HNF4"
                go_part = parts[1].strip()       # "GO:DNA binding ; GO:0003677"

                # Extract InterPro ID
                match = re.search(r'InterPro:(IPR\d+)', ipr_part)
                if not match:
                    continue
                ipr_id = match.group(1)

                # Extract all GO IDs from the line using regex
                go_ids = re.findall(r'GO:\d{7}', go_part)

                # Add GO terms to dictionary
                if go_ids:
                    interpro2go.setdefault(ipr_id, set()).update(go_ids)


    def is_interpro_enzyme(ipr_id, interpro2go, enzyme_go_terms):
        go_terms = interpro2go.get(ipr_id, set())
        return any(go in enzyme_go_terms for go in go_terms)
    
    df['is_enzyme'] = df['InterPro ID'].apply(
        lambda ipr: is_interpro_enzyme(ipr, interpro2go, descendants)
    )

    # now we would like to list the domains which lead to identification as an enzyme, for those genes that are enzymes

    # remove duplicates in the DataFrame
    df = df.drop_duplicates()


    # for each gene, if there are rows where is_enzyme is True, we will concatenate the InterPro IDs and descriptions
    # and keep only the rows where is_enzyme is True
    # for genes where is_enzyme is False that column will be empty

    # thanks genAI for this code which I don't quite understand....
    result_rows = []    
    for name, group in df.groupby('Gene stable ID'):
        has_true = group['is_enzyme'].any()
        has_false = (~group['is_enzyme']).any()

        
        if has_true and has_false:
            # Keep only True rows for concatenation
            filtered_group = group[group['is_enzyme'] == True]
        elif has_true:
            # All rows are True
            filtered_group = group
        else:
            # All rows are False - empty concatenation
            filtered_group = group
        
        # use first row of filtered_group for as template
        result_row = filtered_group.iloc[0].copy()

        if has_true:
            # Create concatenated values and join from true rows
            true_rows = group[group['is_enzyme'] == True]
            concat_values = true_rows['InterPro ID'] + '_' + true_rows['InterPro description']
            result_row['ipID_desc'] = ', '.join(concat_values)
        else:
            result_row['ipID_desc'] = ''
        
        result_rows.append(result_row)

    # Create final DataFrame from result rows
    out = pd.DataFrame(result_rows).reset_index(drop=True)

    out.drop(columns=['InterPro ID', 'InterPro description'], inplace=True)

    return out


