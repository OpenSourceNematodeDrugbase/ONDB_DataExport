import requests
import pandas as pd
from io import StringIO

def fetch_wbp_biomart_using_xml(xml_query):
	"""
	Fetches data from the Wormbase Parasite BioMart using an XML query.
	
	Parameters:
	xml_query (str): The XML query string to be sent to the BioMart service.
	
	Returns:
	pd.DataFrame: A DataFrame containing the results of the query.
	"""
	# Wormbase Parasite BioMart URL
	biomart_url = "https://parasite.wormbase.org/biomart/martservice?query="
	
	# Remove newlines in the xml_query string
	xml_query = xml_query.replace('\n', '').replace('\r', '')
	
	# Change header = "0" to header = "1" to get the column names in the first row of the biomart output
	import re
	xml_query = re.sub(r'header\s*=\s*"0"', 'header="1"', xml_query)

	query_url = biomart_url + xml_query

	# Send the request to the BioMart service
	response = requests.get(query_url)

	if response.status_code == 200:
		return pd.read_csv(StringIO(response.text), sep="\t")
	else:
		response.raise_for_status()

