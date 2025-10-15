"""
populateGeneList.py (v0.3)
Enhanced version with single-species processing, header validation, timeouts, and robust error handling.
This module takes a single parasite species and returns a clean gene table with exactly one row per gene, 
selecting the "main" version (isoform) of each gene.
"""

from __future__ import annotations

import time
from typing import Dict, Set
import pandas as pd
import polars as pl

# Import your existing querybiomart function
from queryWbpBiomart import fetch_biomart_with_validation

# ===== CONFIGURATION =====
REQUEST_TIMEOUT = 120  # seconds
RATE_LIMIT_DELAY = 0.5  # seconds between API calls
MAX_RETRIES = 3

# Expected column names from BioMart (with headers) - MATCHING WHAT BIOMART ACTUALLY RETURNS
EXPECTED_COLUMNS = {
    'wbps_gene_id',
    'wbps_transcript_id', 
    'genome_name',
    'description',
    'transcript_biotype'
}

CDS_EXPECTED_COLUMNS = {
    'wbps_gene_id',
    'wbps_transcript_id', 
    'CDS start (within cDNA)',
    'CDS end (within cDNA)'
}

# Species-specific minimum gene counts
SPECIES_MIN_GENES = {
    'caenorhabditis_elegans': 15000,
    'schistosoma_mansoni': 10000,
    'default': 5000  # Fallback for unknown species
}

# ===== XML TEMPLATES (WITH HEADERS) =====
GENES_XML_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="parasite_mart" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">
    <Dataset name="wbps_gene" interface="default">
        <Filter name="species_id_1010" value="{species_code}"/>
        <Filter name="biotype" value="protein_coding"/>
        <Attribute name="genome_name" />
        <Attribute name="wbps_gene_id" />
        <Attribute name="wbps_transcript_id" />
        <Attribute name="transcript_biotype" />
        <Attribute name="description" />
    </Dataset>
</Query>"""

CDS_XML_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="parasite_mart" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">
    <Dataset name="wbps_gene" interface="default">
        <Filter name="species_id_1010" value="{species_code}"/>
        <Filter name="biotype" value="protein_coding"/>
        <Attribute name="wbps_gene_id" />
        <Attribute name="wbps_transcript_id" />
        <Attribute name="cds_start" />
        <Attribute name="cds_end" />
    </Dataset>
</Query>"""

# ===== CORE FUNCTIONS =====
def validate_biomart_columns(df: pl.DataFrame, expected_columns: Set[str]) -> None:
    """
    Validate that BioMart returned the expected columns.
    Provides clear error messages if headers change.
    """
    if df.is_empty():
        raise ValueError("BioMart returned empty response")
    
    actual_columns = set(df.columns)
    missing_columns = expected_columns - actual_columns
    
    if missing_columns:
        available_columns = ", ".join(sorted(actual_columns))
        raise ValueError(
            f"BioMart response missing expected columns: {sorted(missing_columns)}\n"
            f"Available columns: {available_columns}\n"
            f"This may indicate a BioMart schema change."
        )

def fetch_with_retry(xml_query: str, expected_columns: Set[str], timeout: int = REQUEST_TIMEOUT) -> pl.DataFrame:
    """
    Fetch data from BioMart with retry logic, rate limiting, and header validation.
    """
    for attempt in range(MAX_RETRIES):
        try:
            # Rate limiting between attempts
            if attempt > 0:
                time.sleep(RATE_LIMIT_DELAY * (2 ** attempt))  # Exponential backoff
            
            print(f"Attempt {attempt + 1}/{MAX_RETRIES} for BioMart query...")
            df = fetch_biomart_with_validation(xml_query, expected_columns, timeout)
            
            # Debug output
            print("=== DEBUG COLUMNS ===")
            print("Expected:", expected_columns)
            print("Actual:", set(df.columns))
            print("=====================")
            
            # Validate headers immediately
            validate_biomart_columns(df, expected_columns)
            
            return df
            
        except Exception as e:
            if attempt == MAX_RETRIES - 1:
                raise RuntimeError(f"BioMart query failed after {MAX_RETRIES} attempts: {e}") from e
            print(f"Attempt {attempt + 1} failed, retrying...: {e}")
    
    raise RuntimeError("Unexpected error in fetch_with_retry")

def select_principal_isoform(genes_df: pl.DataFrame, cds_df: pl.DataFrame) -> pl.DataFrame:
    """
    Select principal isoform for each gene based on longest CDS.
    Improved labeling and validation.
    """
    # Calculate CDS length for each transcript
    cds_lengths = cds_df.with_columns(
        (pl.col('CDS end (within cDNA)') - pl.col('CDS start (within cDNA)')).alias('cds_length')
    ).select(['wbps_transcript_id', 'cds_length'])
    
    # Join CDS lengths to genes
    merged = genes_df.join(cds_lengths, on='wbps_transcript_id', how='left')
    
    # Select principal isoform (longest CDS)
    principal_isoforms = (
        merged
        .with_columns(
            # Improved labeling
            pl.when(pl.col('cds_length').is_null())
            .then(pl.lit('non_coding'))
            .when(pl.col('cds_length') == 0)
            .then(pl.lit('non_coding'))
            .otherwise(pl.lit('protein_coding'))
            .alias('isoform_type')
        )
        .sort(['cds_length'], descending=True, nulls_last=True)
        .group_by('wbps_gene_id')
        .agg([
            pl.first('wbps_transcript_id').alias('principal_transcript_id'),
            pl.first('cds_length').alias('principal_cds_length'),
            pl.first('isoform_type').alias('principal_isoform_type'),
            pl.first('description').alias('gene_description'),
            pl.first('genome_name').alias('genome_name'),
            pl.first('transcript_biotype').alias('transcript_biotype')
        ])
    )
    
    return principal_isoforms

def validate_gene_list(df: pl.DataFrame, species_code: str) -> None:
    """
    Comprehensive validation of the final gene list.
    """
    # 1. Check for duplicate Gene IDs (CRITICAL)
    duplicate_count = df.filter(pl.col('wbps_gene_id').is_duplicated()).height
    if duplicate_count > 0:
        raise ValueError(f"Found {duplicate_count} duplicate Gene IDs after principal isoform selection")
    
    # 2. Species-appropriate size validation
    min_genes = SPECIES_MIN_GENES.get(species_code, SPECIES_MIN_GENES['default'])
    if df.height < min_genes:
        raise ValueError(
            f"Only {df.height} genes retrieved for {species_code}. "
            f"Expected at least {min_genes}. Check species code and data availability."
        )
    
    # 3. Validate principal isoform selection worked
    null_principal = df.filter(pl.col('principal_transcript_id').is_null()).height
    if null_principal > 0:
        raise ValueError(f"{null_principal} genes missing principal transcript assignment")
    
    # 4. Check for reasonable CDS lengths in protein-coding genes
    protein_coding = df.filter(pl.col('principal_isoform_type') == 'protein_coding')
    if protein_coding.height > 0:
        short_cds = protein_coding.filter(pl.col('principal_cds_length') < 100).height
        if short_cds > protein_coding.height * 0.1:  # More than 10% have very short CDS
            print(f"Warning: {short_cds} protein-coding genes have very short CDS (<100 bp)")

# ===== MAIN FUNCTION =====
def retrieveGeneListFromWbpBiomart(species_code: str) -> pd.DataFrame:
    """
    Fetch protein-coding genes for a single species and select principal isoforms.
    
    Parameters
    ----------
    species_code : str
        Single WBPS species identifier (e.g., 'wubancprjna275548')
        Do NOT pass comma-separated lists.
        
    Returns
    -------
    pandas.DataFrame
        One row per gene with principal isoform selected by longest CDS.
        Returns pandas DataFrame for compatibility with existing pipelines.
    """
    
    # 1. Input validation - SINGLE species only
    if not isinstance(species_code, str) or not species_code.strip():
        raise ValueError("species_code must be a non-empty string")
    
    if "," in species_code:
        raise ValueError(
            "This function processes ONE species at a time. "
            "Call it multiple times for multiple species."
        )
    
    species_code = species_code.strip()
    print(f"Processing species: {species_code}")
    
    try:
        # 2. Fetch gene data with headers and validation
        genes_xml = GENES_XML_TEMPLATE.format(species_code=species_code)
        genes_df = fetch_with_retry(genes_xml, EXPECTED_COLUMNS)
        
        # 3. Fetch CDS data with headers and validation
        time.sleep(RATE_LIMIT_DELAY)  # Be polite to the server
        cds_xml = CDS_XML_TEMPLATE.format(species_code=species_code)
        cds_df = fetch_with_retry(cds_xml, CDS_EXPECTED_COLUMNS)
        
        # 4. Select principal isoforms (all operations in Polars)
        principal_isoforms = select_principal_isoform(genes_df, cds_df)
        
        # 5. Comprehensive validation
        validate_gene_list(principal_isoforms, species_code)
        
        # 6. Convert to pandas ONLY at the end
        result_df = principal_isoforms.to_pandas()
        
        print(f"Successfully processed {len(result_df)} genes for {species_code}")
        return result_df
        
    except Exception as e:
        raise RuntimeError(f"Failed to process species '{species_code}': {e}") from e


# ===== USAGE EXAMPLE =====
if __name__ == "__main__":
    # Example usage (your main script will call this)
    try:
        result = retrieveGeneListFromWbpBiomart("caenorhabditis_elegans")
        print(f"Retrieved {len(result)} genes")
        print(result.head())
    except Exception as e:
        print(f"Error: {e}")