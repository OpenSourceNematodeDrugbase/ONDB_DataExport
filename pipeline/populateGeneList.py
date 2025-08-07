import pandas as pd

from queryWbpBiomart import *

def retrieveGeneListFromWbpBiomart():

    # define list of genomes to query

    # for now, we will just use the most recent Wuchereria bancrofti genome
    genomes = "wubancprjna275548" 

    # define the XML query to fetch gene IDs and descriptions

    xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                
        <Dataset name = "wbps_gene" interface = "default" >
            <Filter name = "species_id_1010" value = "{genomes}"/>
            <Attribute name = "wbps_gene_id" />
            <Attribute name = "transcript_biotype" />
            <Attribute name = "description" />
        </Dataset>
    </Query>
    """

    # Fetch the data using the XML query
    df = fetch_wbp_biomart_using_xml(xml_query)

    # Tests

    # There should be no dupliacates in the gene IDs

    duplicates = df['Gene stable ID'].duplicated().any()
    assert not duplicates, "Duplicate values found in 'id' column"

    # There should be at least 10000 genes retrieved
    assert len(df) >= 10000, "Less than 10000 genes retrieved"

    return df










