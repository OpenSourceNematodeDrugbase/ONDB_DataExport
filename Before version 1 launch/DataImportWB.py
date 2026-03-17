import requests
import pandas as pd
from io import StringIO

# Wormbase Parasite BioMart URL
biomart_url = "https://parasite.wormbase.org/biomart/martservice?query="

xml_query = """
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "wbps_gene" interface = "default" >
		<Filter name = "species_id_1010" value = "wubancprjna275548"/>
		<Attribute name = "production_name_1010" />
		<Attribute name = "wbps_gene_id" />
		<Attribute name = "hsapiens_gene" />
		<Attribute name = "hsapiens_gene_name" />
		<Attribute name = "hsapiens_homolog_ensembl_peptide" />
		<Attribute name = "hsapiens_orthology_type" />
		<Attribute name = "hsapiens_homolog_perc_id" />
		<Attribute name = "hsapiens_homolog_perc_id_r1" />
	</Dataset>
</Query>
"""

# remove newlines in the xml_query string
xml_query = xml_query.replace('\n', '').replace('\r', '')

# change header = "0" to header = "1" to get the column names in the first row of the biomart output
import re

xml_query = re.sub(r'header\s*=\s*"0"', 'header="1"', xml_query)

import requests

response = requests.get(biomart_url + xml_query)

df = pd.read_csv(StringIO(response.text), sep="\t")

print(df)