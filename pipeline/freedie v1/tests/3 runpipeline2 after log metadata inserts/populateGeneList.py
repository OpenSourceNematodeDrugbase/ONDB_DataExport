import pandas as pd
import polars as pl

from queryWbpBiomart import *

def retrieveGeneListFromWbpBiomart(genomes):

    # define the XML query to fetch gene IDs and descriptions

    xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                
        <Dataset name = "wbps_gene" interface = "default" >
            <Filter name = "species_id_1010" value = "{genomes}"/>
            <Filter name = "biotype" value = "protein_coding"/>
            <Attribute name = "display_name_1010" />
            <Attribute name = "wbps_gene_id" />
            <Attribute name = "wbps_transcript_id" />
            <Attribute name = "transcript_biotype" />
            <Attribute name = "description" />
        </Dataset>
    </Query>
    """

    # Fetch the data using the XML query
    df = fetch_wbp_biomart_using_xml_polars(xml_query)

    # now figure out which transcript to treat as the principal isoform of the gene
    # by fetching exon data for the transcripts

    # one thing to note is that this query is extremely poorly performant on WBP biomart if you ask for multiple genomes
    # time for wb genome = 3 s
    # time to tm genome =  3 s
    # time for both genomes = 380 s
    # hence I'm going to break up this query....

    import time

    genomes_list = genomes.split(",")
    all_exons_data = pl.DataFrame()
    for genome in genomes_list:
        xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                <Dataset name = "wbps_gene" interface = "default" >
                    <Filter name = "species_id_1010" value = "{genome}"/>
                    <Filter name = "biotype" value = "protein_coding"/>
                    <Attribute name = "display_name_1010" />
                    <Attribute name = "wbps_gene_id" />
                    <Attribute name = "wbps_transcript_id" />
                    <Attribute name = "cds_start" />
                    <Attribute name = "cds_end" />
                </Dataset>
            </Query>
        """    
        exons_data = fetch_wbp_biomart_using_xml_polars(xml_query)
        all_exons_data = pl.concat([all_exons_data, exons_data], how="vertical")

    all_exons_data2 = all_exons_data.group_by('Transcript stable ID').agg(pl.col('CDS end (within cDNA)').max().alias('CDS_max'))

    df2 = df.join(all_exons_data2, left_on='Transcript stable ID', right_on='Transcript stable ID', how='left') 

    # now we can find the principal isoform of each gene by taking the transcript with the longest CDS
    # nulls_last=True because we don't want to choose a non-coding isoform
    result = df2.sort('CDS_max', descending=True, nulls_last=True).group_by('Gene stable ID', maintain_order=True).first()                     

    result_ordered = result.sort('Genome name', 'Gene stable ID')


    df = result_ordered.to_pandas()

    # Tests

    # There should be no duplicates in the gene IDs

    duplicates = df['Gene stable ID'].duplicated().any()
    assert not duplicates, "Duplicate values found in 'id' column"

    # There should be at least 10000 genes retrieved
    assert len(df) >= 10000, "Less than 10000 genes retrieved"

    return df










