import pandas as pd
from queryWbpBiomart import *

def queryWbpHumanOrthologues(genomes):

    xml_query = f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >		
	<Dataset name = "wbps_gene" interface = "default" >
		<Filter name = "species_id_1010" value = "{genomes}"/>
		<Attribute name = "wbps_gene_id" />
		<Attribute name = "hsapiens_gene" />
		<Attribute name = "hsapiens_gene_name" />
		<Attribute name = "hsapiens_orthology_type" />
		<Attribute name = "hsapiens_homolog_perc_id" />
		<Attribute name = "hsapiens_homolog_perc_id_r1" />
	</Dataset>
    </Query>
    """

    df = fetch_wbp_biomart_using_xml(xml_query)

    # now we have an issue here that there can be multiple human orthologues for a single parasite gene
    # so we will look for duplicated rows (same gene ID) and keep the rows with the highest percentage identity score
    # as the higher the percentage identity the worse the situation for our purpose (wanting to avoid similarity to human genes)
    df['% identity'] = df['% identity'].fillna(0)
    df_filtered = df.loc[df.groupby('Gene stable ID')['% identity'].idxmax()]


    # Apply criterion "lacks WBP human orthologue"
    df_filtered['lacks_WBP_human_orthologue'] = pd.isna(df_filtered['Human gene stable ID'])

    # Apply criterion "best WBP human orthologue < 40% identity"
    df_filtered['best_WBP_human_orthologue_lt_40pct_identity'] = df_filtered['% identity'] < 40



    # Tests
    # There should be no dupliacates in the gene IDs
    duplicates = df_filtered['Gene stable ID'].duplicated().any()
    assert not duplicates, "Duplicate values found in 'id' column"
    # There should be at least 10000 genes retrieved

    assert len(df_filtered) >= 10000, "Less than 10000 genes retrieved"



    return df_filtered


