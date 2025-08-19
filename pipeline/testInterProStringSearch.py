from numpy import sort
import polars as pl
from queryWbpBiomart import fetch_wbp_biomart_using_xml_polars


def testInterProStringSearch(genomes, search_string, test_name):
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

    # concatenate interpro_id and interpro_description
    df = df.with_columns((pl.col("InterPro ID") + " - " + pl.col("InterPro description")).alias("interpro_annotation"))

    q = df.with_columns(
        # can't use str.contains.any as it does not support regular expressions
        pl.when(pl.col("interpro_annotation").str.contains(search_string))
        .then(pl.col("interpro_annotation"))
        .otherwise(pl.lit(None))
        .alias('matches')
    )
    
    match = q.group_by("Gene stable ID").agg(pl.col("matches").str.join(", ").alias(test_name+'_evidence')).sort('Gene stable ID')

    output = match.with_columns(
        pl.when(pl.col(test_name+'_evidence') == '')
        .then(pl.lit(False))
        .otherwise(pl.lit(True))
        .alias(test_name)
    )    

    output = output.with_columns(
       pl.when(pl.col(test_name) == False)
         .then(pl.lit("Encodes protein that lacks InterPro GPCR domains"))
         .otherwise('Encodes protein with InterPro domain(s): ' + pl.col(test_name+'_evidence'))
         .alias(test_name+'_evidence')
    )

    return(output)



