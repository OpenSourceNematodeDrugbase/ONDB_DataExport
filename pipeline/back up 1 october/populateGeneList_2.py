"""
populateGeneList.py (v0.2)
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
HEADER_AUDIT = {}



EXPECTED_COLUMNS = {
    'Genome project',
    'Gene stable ID',
    'Transcript stable ID', 
    'Protein stable ID',
    'Gene name',
    'Gene description',
    'Gene biotype',
    'Chromosome/scaffold name',
    'Gene start (bp)', 
    'Gene end (bp)',
    'Strand',
    'Transcript start (bp)',
    'Transcript end (bp)',
    '% GC content',
    'Transcript count'
}

CDS_EXPECTED_COLUMNS = {
    'Genome project',
    'Gene stable ID',
    'Transcript stable ID', 
    'Protein stable ID',
    'CDS start (within cDNA)',
    'CDS end (within cDNA)',
    'Transcript start (bp)',
    'Transcript end (bp)',
    'Chromosome/scaffold name',
    'Gene start (bp)',
    'Gene end (bp)',
    'Strand',
    'Gene biotype',
    'Gene description',
    '% GC content',
    'Transcript count'
}


## WORKING 28/09/2025 - CHECK POINT
# EXPECTED_COLUMNS = {
#     'Gene stable ID',
#     'Transcript stable ID', 
#     'Gene name',
#     'Gene description',     
#     'Gene biotype'    
# }

# CDS_EXPECTED_COLUMNS = {
#     'Gene stable ID',
#     'Transcript stable ID', 
#     'CDS start (within cDNA)',    
#     'CDS end (within cDNA)'     
# }



## CHECK WHAT THIS DOES?

# Species-specific minimum gene counts
SPECIES_MIN_GENES = {
    'caenorhabditis_elegans': 15000,
    'schistosoma_mansoni': 10000,
    'default': 5000  # Fallback for unknown species
}

# === Phase 2: normalization state ===
# Will be written into run.json at the end of the run (from run_v01.py)
NORMALIZE_WARNINGS = {
    "missing_columns": {},   # species -> [cols]
    "type_coercions": {}     # species -> {col: count}
}

# Final canonical export order (stable API for downstream)
CANONICAL_EXPORT_ORDER = [
    "wbps_gene_id",
    "principal_transcript_id",
    "principal_cds_length",
    "principal_isoform_type",
    "genome_name",
    "gene_description",
    "transcript_biotype",
    "species",
]

## NEW VERSION

# CANONICAL_EXPORT_ORDER = [
#     "wbps_gene_id",
#     "principal_transcript_id", 
#     "principal_cds_length",
#     "principal_isoform_type",
#     "gene_name",            # 🆕 NEW NAME
#     "gene_description",     # 🆕 KEEPS SAME FINAL NAME
#     "transcript_biotype",   # 🆕 KEEPS SAME FINAL NAME  
#     "species",
# ]


# ===== XML TEMPLATES (WITH HEADERS) =====        <Filter name="gene_biotype" value="protein_coding"/>




## EXPANDING GRADUALLY
GENES_XML_TEMPLATE = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="parasite_mart" formatter="TSV" header="1" uniqueRows="0">
    <Dataset name="wbps_gene" interface="default">
        <Filter name="species_id_1010" value="{species_code}"/>
        <Attribute name="production_name_1010"/>   <!-- Genome project -->
        <Attribute name="wbps_gene_id"/>           <!-- Gene stable ID -->
        <Attribute name="wbps_transcript_id"/>     <!-- Transcript stable ID -->
        <Attribute name="wbps_peptide_id"/>        <!-- Protein stable ID -->
        <Attribute name="external_gene_id"/>       <!-- Gene name -->
        <Attribute name="description"/>            <!-- Gene description -->
        <Attribute name="gene_biotype"/>           <!-- Gene biotype -->
        <Attribute name="chromosome_name"/>        <!-- Chromosome/scaffold name -->
        <Attribute name="start_position"/>         <!-- Gene start (bp) -->
        <Attribute name="end_position"/>           <!-- Gene end (bp) -->
        <Attribute name="strand"/>                 <!-- Strand -->
        <Attribute name="transcript_start"/>       <!-- Transcript start (bp) -->
        <Attribute name="transcript_end"/>         <!-- Transcript end (bp) -->
        <Attribute name="percentage_gc_content"/>  <!-- % GC content -->
        <Attribute name="transcript_count"/>       <!-- Transcript count -->
    </Dataset>
</Query>'''


## CDS GRADUALLY BIGGER TEMPLATE

CDS_XML_TEMPLATE = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="parasite_mart" formatter="TSV" header="1" uniqueRows="0">
    <Dataset name="wbps_gene" interface="default">
        <Filter name="species_id_1010" value="{species_code}"/>
        <Filter name="biotype" value="protein_coding"/>
        <Attribute name="production_name_1010"/>   <!-- Genome project -->
        <Attribute name="wbps_gene_id"/>           <!-- Gene stable ID -->
        <Attribute name="wbps_transcript_id"/>     <!-- Transcript stable ID -->
        <Attribute name="wbps_peptide_id"/>        <!-- Protein stable ID -->
        <Attribute name="cds_start"/>              <!-- CDS start (within cDNA) -->
        <Attribute name="cds_end"/>                <!-- CDS end (within cDNA) -->
        <Attribute name="transcript_start"/>       <!-- Transcript start (bp) -->
        <Attribute name="transcript_end"/>         <!-- Transcript end (bp) -->
        <Attribute name="chromosome_name"/>        <!-- Chromosome/scaffold name -->
        <Attribute name="start_position"/>         <!-- Gene start (bp) -->
        <Attribute name="end_position"/>           <!-- Gene end (bp) -->
        <Attribute name="strand"/>                 <!-- Strand -->
        <Attribute name="gene_biotype"/>           <!-- Gene biotype -->
        <Attribute name="description"/>            <!-- Gene description -->
        <Attribute name="percentage_gc_content"/>  <!-- % GC content -->
        <Attribute name="transcript_count"/>       <!-- Transcript count -->
    </Dataset>
</Query>'''



#### check point 28/09/2025 ###


###WORKING TEMPLATE 28/09/2025

# CDS_XML_TEMPLATE = '''<?xml version="1.0" encoding="UTF-8"?>
# <!DOCTYPE Query>
# <Query virtualSchemaName="parasite_mart" formatter="TSV" header="1" uniqueRows="0">
#     <Dataset name="wbps_gene" interface="default">
#         <Filter name="species_id_1010" value="{species_code}"/>
#         <Filter name="biotype" value="protein_coding"/>
#         <Attribute name="wbps_gene_id"/>
#         <Attribute name="wbps_transcript_id"/>
#         <Attribute name="cds_start"/>
#         <Attribute name="cds_end"/>
#     </Dataset>
# </Query>'''


###WORKING TEMPLATE 28/09/2025

# GENES_XML_TEMPLATE = '''<?xml version="1.0" encoding="UTF-8"?>
# <!DOCTYPE Query>
# <Query virtualSchemaName="parasite_mart" formatter="TSV" header="1" uniqueRows="0">
#     <Dataset name="wbps_gene" interface="default">
#         <Filter name="species_id_1010" value="{species_code}"/>
#         <Attribute name="wbps_gene_id"/>           <!-- Gene stable ID -->
#         <Attribute name="wbps_transcript_id"/>     <!-- Transcript stable ID -->
#         <Attribute name="external_gene_id"/>       <!-- Gene name -->
#         <Attribute name="description"/>            <!-- Gene description -->
#         <Attribute name="gene_biotype"/>           <!-- Gene biotype -->
#     </Dataset>
# </Query>
# '''


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

# def select_principal_isoform(genes_df: pl.DataFrame, cds_df: pl.DataFrame) -> pl.DataFrame:
#     """
#     Select principal isoform for each gene based on longest CDS.
#     Improved labeling and validation.
#     """
#     # Calculate CDS length for each transcript
#     cds_lengths = cds_df.with_columns(
#         (pl.col('CDS end (within cDNA)') - pl.col('CDS start (within cDNA)')).alias('cds_length')
#     ).select(['Transcript stable ID', 'cds_length'])
    
#     # Join CDS lengths to genes
#     merged = genes_df.join(cds_lengths, on='Transcript stable ID', how='left')
    
#     # Select principal isoform (longest CDS)
#     principal_isoforms = (
#         merged
#         .with_columns(
#             # Improved labeling
#             pl.when(pl.col('cds_length').is_null())
#             .then(pl.lit('non_coding'))
#             .when(pl.col('cds_length') == 0)
#             .then(pl.lit('non_coding'))
#             .otherwise(pl.lit('protein_coding'))
#             .alias('isoform_type')
#         )
#    #     .sort(['cds_length'], descending=True, nulls_last=True) Non deterministic ISSSUE
#         .sort(['cds_length', 'Gene stable ID', 'Transcript stable ID'], descending=[True, False, False], nulls_last=True)
       
#         # In select_principal_isoform(), right after the sort:
#         .group_by('Gene stable ID')

# ## ORIFINAL NAMES ####

#             .agg([
#         pl.first('Transcript stable ID').alias('principal_transcript_id'),
#         pl.first('cds_length').alias('principal_cds_length'),
#         pl.first('isoform_type').alias('principal_isoform_type'),
#         pl.first('Gene description').alias('gene_description'),           # ← ORIGINAL NAME
#         pl.first('Gene name').alias('genome_name'),                       # ← ORIGINAL NAME
#         pl.first('Gene biotype').alias('transcript_biotype')              # ← ORIGINAL NAME
#     ])
#     )
    
#         # NEW VERSION BUT NEEDS ADJUSTMENT
#         # .agg([
#         #     pl.first('Transcript stable ID').alias('principal_transcript_id'),
#         #     pl.first('cds_length').alias('principal_cds_length'),
#         #     pl.first('isoform_type').alias('principal_isoform_type'),
#         #     pl.first('Gene description').alias('principal_transcript_description'),  # 🆕 NEW NAME
#         #     pl.first('Gene name').alias('gene_name'),                                # 🆕 NEW NAME
#         #     pl.first('Gene biotype').alias('principal_transcript_biotype')           # 🆕 NEW NAME
#         # ])
#     return principal_isoforms




def select_principal_isoform(genes_df: pl.DataFrame, cds_df: pl.DataFrame) -> pl.DataFrame:
    """
    Select principal isoform for each gene based on longest CDS.
    Improved labeling and validation.
    """
    # Calculate CDS length for each transcript
    cds_lengths = cds_df.with_columns(
        (pl.col('CDS end (within cDNA)') - pl.col('CDS start (within cDNA)')).alias('cds_length')
    ).select(['Transcript stable ID', 'cds_length'])
    
    # Join CDS lengths to genes
    merged = genes_df.join(cds_lengths, on='Transcript stable ID', how='left')
    
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
        .sort(['cds_length', 'Gene stable ID', 'Transcript stable ID'], descending=[True, False, False], nulls_last=True)
        .group_by('Gene stable ID')
        .agg([
            pl.first('Transcript stable ID').alias('principal_transcript_id'),
            pl.first('cds_length').alias('principal_cds_length'),
            pl.first('isoform_type').alias('principal_isoform_type'),
            pl.first('Gene description').alias('gene_description'),
            pl.first('Gene name').alias('genome_name'),
            pl.first('Gene biotype').alias('transcript_biotype')
        ])
    )
    
    return principal_isoforms










def validate_gene_list(df: pl.DataFrame, species_code: str) -> None:
    """
    Comprehensive validation of the final gene list.
    """
    # 1. Check for duplicate Gene IDs (CRITICAL)
    duplicate_count = df.filter(pl.col('Gene stable ID').is_duplicated()).height
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



def _normalize_per_species(principal_df: pl.DataFrame, species_code: str) -> pl.DataFrame:
    """
    Rename friendly BioMart labels to canonical names, ensure required columns exist,
    coerce basic types, add `species`, and record any warnings.
    """
    # 1) Rename friendly headers -> canonical

    rename_map = {
    "Gene stable ID": "wbps_gene_id",
    "Transcript stable ID": "wbps_transcript_id",  # include if present
    "Gene name": "genome_name",                 # <- new (for display_name_1010)
    "Genome name": "genome_name",                  # keep, in case some datasets use it / assuming display name referes to genomes
    "Gene description": "gene_description",
    "Transcript biotype": "transcript_biotype",
    "gene_name": "gene_description",               # legacy internal alias, optional
    }

    # rename_map = {
    #     "Gene stable ID": "wbps_gene_id",
    #     "gene_name": "gene_name",                               # 🆕 NEW MAPPING
    #     "principal_transcript_description": "gene_description", # 🆕 NEW MAPPING  
    #     "principal_transcript_biotype": "transcript_biotype",   # 🆕 NEW MAPPING
    # }


    df = principal_df.rename({k: v for k, v in rename_map.items() if k in principal_df.columns})

    # 2) Add species column
    df = df.with_columns(pl.lit(species_code).alias("species"))

    # 3) Ensure all required canonical columns exist; warn & fill with nulls if missing
    missing = [col for col in CANONICAL_EXPORT_ORDER if col not in df.columns]
    if missing:
        df = df.with_columns([pl.lit(None).alias(col) for col in missing])
        NORMALIZE_WARNINGS["missing_columns"][species_code] = missing

    # 4) Basic type coercions (soft)
      # - IDs & text → string
    str_cols = [
        "wbps_gene_id",
        "principal_transcript_id",
        "genome_name",
        "gene_description",
        "transcript_biotype",
        "species",
    ]


# NEW 
    # str_cols = [
    #     "wbps_gene_id",
    #     "principal_transcript_id",
    #     "gene_name",        # 🆕 NEW NAME
    #     "gene_description",
    #     "transcript_biotype",
    #     "species", 
    # ]


    df = df.with_columns([pl.col(c).cast(pl.Utf8) for c in str_cols if c in df.columns])

    #    - Numeric → Int64 (nullable). Count nulls created by coercion (if any).
    if "principal_cds_length" in df.columns:
        before_non_null = int(df.select(pl.col("principal_cds_length").is_not_null().sum()).item())
        df = df.with_columns(pl.col("principal_cds_length").cast(pl.Int64, strict=False))
        after_non_null = int(df.select(pl.col("principal_cds_length").is_not_null().sum()).item())
        coerced_to_null = max(0, before_non_null - after_non_null)
        if coerced_to_null > 0:
            NORMALIZE_WARNINGS["type_coercions"].setdefault(species_code, {})["principal_cds_length"] = coerced_to_null

        # Treat negatives as invalid → null (rare but safer)
        neg_count = int(df.filter(pl.col("principal_cds_length") < 0).height)
        if neg_count > 0:
            df = df.with_columns(
                pl.when(pl.col("principal_cds_length") < 0)
                  .then(pl.lit(None))
                  .otherwise(pl.col("principal_cds_length"))
                  .alias("principal_cds_length")
            )
            NORMALIZE_WARNINGS["type_coercions"].setdefault(species_code, {})["principal_cds_length_negatives"] = neg_count

    # 5) Enforce final column order for downstream
    df = df.select(CANONICAL_EXPORT_ORDER)

    return df


#fetch with retry replaced with this one

def fetch_with_retry(xml_query: str, expected_columns: Set[str], timeout: int = REQUEST_TIMEOUT) -> pl.DataFrame:
    last_err = None
    for attempt in range(MAX_RETRIES):
        try:
            if attempt > 0:
                time.sleep(RATE_LIMIT_DELAY * (2 ** attempt))  # backoff
            df = fetch_biomart_with_validation(xml_query, expected_columns=expected_columns, timeout=timeout)
            print("[headers]", list(df.columns))  # debug/audit
            validate_biomart_columns(df, expected_columns)
            return df
        except Exception as e:
            last_err = e
            if attempt == MAX_RETRIES - 1:
                break
            print(f"Attempt {attempt + 1} failed, retrying...: {e}")
    raise RuntimeError(f"BioMart query failed after {MAX_RETRIES} attempts: {last_err}")


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
        

        print(f"[debug] Testing genes XML for species: {species_code}")
        print(f"[debug] First 500 chars of XML:\n{genes_xml[:200]}...")



        genes_df = fetch_with_retry(genes_xml, EXPECTED_COLUMNS)
        
        # 3. Fetch CDS data with headers and validation
        time.sleep(RATE_LIMIT_DELAY)  # Be polite to the server
        cds_xml = CDS_XML_TEMPLATE.format(species_code=species_code)
        
    

        cds_df = fetch_with_retry(cds_xml, CDS_EXPECTED_COLUMNS)
        print(f"[debug] ACTUAL CDS headers: {cds_df.columns}")  # debugging what actuall names are present



        # record the raw headers BioMart returned (friendly labels; header=1)
        HEADER_AUDIT[species_code] = {
            "genes": list(genes_df.columns),
            "cds": list(cds_df.columns),
        }

        # 4. Select principal isoforms (all operations in Polars)
        principal_isoforms = select_principal_isoform(genes_df, cds_df)
        
        # 5. Comprehensive validation
        validate_gene_list(principal_isoforms, species_code)
        # 5.5. Normalize to canonical export schema (stable names)
        canonical_df = _normalize_per_species(principal_isoforms, species_code)



        # # 6. Convert to pandas ONLY at the end
        # result_df = principal_isoforms.to_pandas()
        
        # print(f"Successfully processed {len(result_df)} genes for {species_code}")
        # return result_df

        # 6. Convert to pandas ONLY at the end
        result_df = canonical_df.to_pandas()

        print(f"Successfully processed {len(result_df)} genes for {species_code}")
        return result_df

        
    except Exception as e:
        raise RuntimeError(f"Failed to process species '{species_code}': {e}") from e



# # ===== USAGE EXAMPLE =====
# if __name__ == "__main__":
#     # Example usage (your main script will call this)
#     try:
#         result = retrieveGeneListFromWbpBiomart("caenorhabditis_elegans")
#         print(f"Retrieved {len(result)} genes")
#     except Exception as e:
#         print(f"Error: {e}")