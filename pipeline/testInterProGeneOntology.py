import polars as pl
import pandas as pd
from queryWbpBiomart import *


def testInterProGeneOntology(genomes, go_term, test_name):
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
    df = fetch_wbp_biomart_using_xml_polars(xml_query)

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

    # Find all GO terms descending from the GO term of interest e.g. "GO:0003824" = catalytic activity
    root_go = go_term
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




    def is_go_term_matching(ipr_id, interpro2go, search_go_terms):
        # Check if any GO terms associated with the InterPro ID match the search terms
        go_terms = interpro2go.get(ipr_id, set())
        return any(go in search_go_terms for go in go_terms)
    
    #df[test_name] = df['InterPro ID'].apply(
    #    lambda ipr: is_go_term_matching(ipr, interpro2go, descendants)
    #)

    df = df.with_columns(
        pl.col('InterPro ID')
        .map_elements(lambda ipr: is_go_term_matching(ipr, interpro2go, descendants), return_dtype=pl.Boolean)
        .alias(test_name)
    )

    # now we would like to list the domains which lead to identification as X, for those genes that are X

    # remove duplicates in the DataFrame
    df = df.unique()


    # for each gene, if there are rows where test is True, we will concatenate the InterPro IDs and descriptions
    # and keep only the rows where test is True
    # for genes where test is False that column will be empty

    # concatenate interpro_id and interpro_description
    df = df.with_columns((pl.col("InterPro ID") + " - " + pl.col("InterPro description")).alias("interpro_annotation"))

    qp = df.with_columns(
        # can't use str.contains.any as it does not support regular expressions
        pl.when(pl.col(test_name) == True)
        .then(pl.col("interpro_annotation"))
        .otherwise(pl.lit(None))
        .alias('matches')
    )

    match = qp.group_by("Gene stable ID").agg(pl.col("matches").str.join(", ").alias(test_name+'_evidence')).sort('Gene stable ID')

    output = match.with_columns(
        pl.when(pl.col(test_name+'_evidence') == '')
        .then(pl.lit(False))
        .otherwise(pl.lit(True))
        .alias(test_name)
    )    

    output = output.with_columns(
       pl.when(pl.col(test_name) == False)
         .then(pl.lit("Does not have InterPro domains related to " + go_term))
         .otherwise('Encodes protein with InterPro domain(s)  that are related to ' + go_term + ': ' + pl.col(test_name+'_evidence') )
         .alias(test_name+'_evidence')
    )

    pandas_output = output.to_pandas()

    return pandas_output


