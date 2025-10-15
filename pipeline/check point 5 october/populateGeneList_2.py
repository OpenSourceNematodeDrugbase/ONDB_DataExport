"""
populateGeneList.py (v0.3)
Fixed version with consistent header=1 processing and deterministic isoform selection.
"""

from __future__ import annotations

import time
from typing import Set
import pandas as pd
import polars as pl

from queryWbpBiomart import fetch_biomart_with_validation

# ===== CONFIGURATION =====
REQUEST_TIMEOUT = 120
RATE_LIMIT_DELAY = 0.5
MAX_RETRIES = 3

# Expected columns match EXACT XML attribute names (header="1")
EXPECTED_COLUMNS = {
    'Gene stable ID',
    'Transcript stable ID', 
    'Genome name',
    'Gene description',     
    'Transcript biotype'    
}

CDS_EXPECTED_COLUMNS = {
    'Gene stable ID',
    'Transcript stable ID', 
    'CDS start (within cDNA)',    
    'CDS end (within cDNA)'     
}





# Species-specific minimum gene counts
SPECIES_MIN_GENES = {
    'caenorhabditis_elegans': 15000,
    'schistosoma_mansoni': 10000,
    'default': 5000
}

# # Final canonical export order not used yet
# CANONICAL_EXPORT_ORDER = [
#     "wbps_gene_id",
#     "principal_transcript_id", 
#     "principal_cds_length",
#     "principal_isoform_type",
#     "genome_name",
#     "gene_description",
#     "transcript_biotype",
#     "species",
# ]

# XML Templates (header="1" - raw attribute names)
GENES_XML_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="parasite_mart" formatter="TSV" header="1" uniqueRows="0">
    <Dataset name="wbps_gene" interface="default">
        <Filter name="species_id_1010" value="{species_code}"/>
        <Filter name="biotype" value="protein_coding"/>
        <Attribute name="display_name_1010"/>
        <Attribute name="wbps_gene_id"/>
        <Attribute name="wbps_transcript_id"/>
        <Attribute name="transcript_biotype"/>
        <Attribute name="description"/>
    </Dataset>
</Query>"""

CDS_XML_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="parasite_mart" formatter="TSV" header="1" uniqueRows="0">
    <Dataset name="wbps_gene" interface="default">
        <Filter name="species_id_1010" value="{species_code}"/>
        <Filter name="biotype" value="protein_coding"/>
        <Attribute name="display_name_1010"/>
        <Attribute name="wbps_gene_id"/>
        <Attribute name="wbps_transcript_id"/>
        <Attribute name="cds_start"/>
        <Attribute name="cds_end"/>
    </Dataset>
</Query>"""

# ===== CORE FUNCTIONS =====

def validate_biomart_columns(df: pl.DataFrame, expected_columns: Set[str]) -> None:
    """Validate that BioMart returned the expected columns."""
    if df.is_empty():
        raise ValueError("BioMart returned empty response")
    
    actual_columns = set(df.columns)
    missing_columns = expected_columns - actual_columns
    
    if missing_columns:
        available_columns = ", ".join(sorted(actual_columns))
        raise ValueError(
            f"BioMart response missing expected columns: {sorted(missing_columns)}\n"
            f"Available columns: {available_columns}"
        )

def select_principal_isoform(genes_df: pl.DataFrame, cds_df: pl.DataFrame) -> pl.DataFrame:
    """
    Select principal isoform for each gene based on longest CDS.
    Uses raw attribute names consistently throughout.
    """
   
    # Calculate CDS length using ACTUAL column names
    cds_lengths = cds_df.with_columns(
        (pl.col('CDS end (within cDNA)') - pl.col('CDS start (within cDNA)')).alias('cds_length')
    ).select(['Transcript stable ID', 'cds_length'])

    # Join CDS lengths to genes using ACTUAL column names
    merged = genes_df.join(cds_lengths, on='Transcript stable ID', how='left')  # ← FIXED

    # Select principal isoform with deterministic sorting
    principal_isoforms = (
        merged
        .with_columns(
            pl.when(pl.col('cds_length').is_null())
            .then(pl.lit('non_coding'))
            .when(pl.col('cds_length') == 0)
            .then(pl.lit('non_coding'))
            .otherwise(pl.lit('protein_coding'))
            .alias('principal_isoform_type')
        )
        # Deterministic sort with ACTUAL column names
        .sort(['cds_length', 'Gene stable ID', 'Transcript stable ID'],  
              descending=[True, False, False], nulls_last=True)
        .group_by('Gene stable ID')  # ← FIXED
        # Create final column names using ACTUAL source columns
        .agg([
            pl.first('Transcript stable ID').alias('principal_transcript_id'),
            pl.first('cds_length').alias('principal_cds_length'),
            pl.first('principal_isoform_type').alias('principal_isoform_type'),
            pl.first('Gene description').alias('gene_description'),  
            pl.first('Genome name').alias('genome_name'),  
            pl.first('Transcript biotype').alias('transcript_biotype')  
        ])
    )
    return principal_isoforms.sort(['Gene stable ID'])  

    

def validate_gene_list(df: pl.DataFrame, species_code: str) -> None:
    """Comprehensive validation using correct column names."""
    # Check for duplicate Gene IDs
    duplicate_count = df.filter(pl.col('Gene stable ID').is_duplicated()).height
    if duplicate_count > 0:
        raise ValueError(f"Found {duplicate_count} duplicate Gene IDs")
    
    # Species-appropriate size validation
    min_genes = SPECIES_MIN_GENES.get(species_code, SPECIES_MIN_GENES['default'])
    if df.height < min_genes:
        raise ValueError(f"Only {df.height} genes retrieved for {species_code}. Expected at least {min_genes}")
    
    # Validate principal isoform selection
    null_principal = df.filter(pl.col('principal_transcript_id').is_null()).height
    if null_principal > 0:
        raise ValueError(f"{null_principal} genes missing principal transcript")


def fetch_with_retry(xml_query: str, expected_columns: Set[str], timeout: int = REQUEST_TIMEOUT) -> pl.DataFrame:
    last_err = None
    for attempt in range(MAX_RETRIES):
        try:
            if attempt > 0:
                time.sleep(RATE_LIMIT_DELAY * (2 ** attempt))
            df = fetch_biomart_with_validation(xml_query, expected_columns=expected_columns, timeout=timeout)
            
            # DEBUG: Show actual vs expected
            print(f"[DEBUG] ACTUAL COLUMNS: {df.columns}")
            print(f"[DEBUG] EXPECTED COLUMNS: {expected_columns}")
            
            validate_biomart_columns(df, expected_columns)
            return df
        except Exception as e:
            last_err = e
            print(f"[DEBUG] ERROR: {e}")  # Show the exact error
            if attempt == MAX_RETRIES - 1:
                break
            print(f"Attempt {attempt + 1} failed, retrying...: {e}")
    raise RuntimeError(f"BioMart query failed after {MAX_RETRIES} attempts: {last_err}")


# ===== MAIN FUNCTION =====
def retrieveGeneListFromWbpBiomart(species_code: str) -> pd.DataFrame:
    """
    Fetch protein-coding genes for a single species and select principal isoforms.
    """
    if not isinstance(species_code, str) or not species_code.strip():
        raise ValueError("species_code must be a non-empty string")
    
    if "," in species_code:
        raise ValueError("This function processes ONE species at a time.")
    
    species_code = species_code.strip()
    print(f"Processing species: {species_code}")
    
    try:
        # Fetch gene data
        genes_xml = GENES_XML_TEMPLATE.format(species_code=species_code)
        genes_df = fetch_with_retry(genes_xml, EXPECTED_COLUMNS)
        print(f"[info] BioMart (gene list) — received columns: {list(genes_df.columns)}")
        print(f"[info] BioMart (gene list) — received rows: {genes_df.height:,}")


        # Fetch CDS data
        time.sleep(RATE_LIMIT_DELAY)
        cds_xml = CDS_XML_TEMPLATE.format(species_code=species_code)
        cds_df = fetch_with_retry(cds_xml, CDS_EXPECTED_COLUMNS)
        print(f"[info] BioMart (CDS) — received columns: {list(cds_df.columns)}")
        print(f"[info] BioMart (CDS) — received rows: {cds_df.height:,}")


        # Select principal isoforms
        principal_isoforms = select_principal_isoform(genes_df, cds_df)
        
        # Add species column and validate
        result_df = principal_isoforms.with_columns(pl.lit(species_code).alias("species"))
        validate_gene_list(result_df, species_code)
        
        # Convert to pandas
        final_df = result_df.to_pandas()
        print(f"Successfully processed {len(final_df)} genes for {species_code}")
        return final_df
        
    except Exception as e:
        raise RuntimeError(f"Failed to process species '{species_code}': {e}") from e