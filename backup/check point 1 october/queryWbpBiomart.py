"""
queryWbpBiomart.py (v0.2)
Enhanced BioMart client with timeouts, retries, and proper error handling.
Aligns with populateGeneList.py framework.
"""

from __future__ import annotations

import re
import time
from urllib.parse import quote
from io import StringIO
from typing import Union

import requests
import pandas as pd
import polars as pl

# Configuration aligned with populateGeneList.py
BIOMART_URL = "https://parasite.wormbase.org/biomart/martservice"
REQUEST_TIMEOUT = 120  # seconds
MAX_RETRIES = 3
RATE_LIMIT_DELAY = 0.5  # seconds

def fetch_wbp_biomart(
    xml_query: str, 
    engine: str = "polars",
    timeout: int = REQUEST_TIMEOUT,
    max_retries: int = MAX_RETRIES
) -> Union[pd.DataFrame, pl.DataFrame]:
    """
    Unified function to fetch data from WormBase ParaSite BioMart.
    
    Parameters
    ----------
    xml_query : str
        The XML query string for BioMart
    engine : str, default "polars"
        Return format: "polars" or "pandas"
    timeout : int, default 120
        Request timeout in seconds
    max_retries : int, default 3
        Maximum number of retry attempts
        
    Returns
    -------
    Union[pd.DataFrame, pl.DataFrame]
        DataFrame in the requested format
    """
    # Validate inputs
    if engine not in ["polars", "pandas"]:
        raise ValueError("engine must be 'polars' or 'pandas'")
    
    if not xml_query or not isinstance(xml_query, str):
        raise ValueError("xml_query must be a non-empty string")
    
   
    
    # Prepare the query URL safely
    query_url = _build_biomart_url(xml_query)
    
    # Fetch with retry logic
    response_text = _fetch_with_retry(query_url, timeout, max_retries)
    
    # Parse based on requested engine
    if engine == "polars":
        return pl.read_csv(StringIO(response_text), separator="\t")
    else:
        return pd.read_csv(StringIO(response_text), sep="\t")






def _build_biomart_url(xml_query: str) -> str:
    """Safely construct the BioMart URL with proper encoding."""
    # URL-encode the XML query to handle special characters
    encoded_query = quote(xml_query)
    return f"{BIOMART_URL}?query={encoded_query}"

def _fetch_with_retry(url: str, timeout: int, max_retries: int) -> str:
    """Fetch data with retry logic and rate limiting."""
    last_exception = None
    
    for attempt in range(max_retries):
        try:
            # Rate limiting between retries
            if attempt > 0:
                time.sleep(RATE_LIMIT_DELAY * (2 ** attempt))
            
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()
            
            # Check if response is not empty
            if not response.text.strip():
                raise ValueError("BioMart returned empty response")
                
            return response.text
            
        except requests.exceptions.Timeout as e:
            last_exception = e
            print(f"Timeout on attempt {attempt + 1}/{max_retries}")
        except requests.exceptions.HTTPError as e:
            last_exception = e
            if e.response.status_code >= 500:  # Server errors can be retried
                print(f"Server error {e.response.status_code} on attempt {attempt + 1}")
            else:
                raise  # Client errors (4xx) shouldn't be retried
        except Exception as e:
            last_exception = e
            print(f"Error on attempt {attempt + 1}/{max_retries}: {e}")
    
    # All retries failed
    raise RuntimeError(f"Failed after {max_retries} attempts: {last_exception}")

# Backward compatibility functions
# def fetch_wbp_biomart_using_xml(xml_query: str, timeout: int = REQUEST_TIMEOUT) -> pd.DataFrame:
#     """Legacy function for pandas output."""
#     return fetch_wbp_biomart(xml_query, engine="pandas", timeout=timeout)

def fetch_wbp_biomart_using_xml_polars(xml_query: str, timeout: int = REQUEST_TIMEOUT) -> pl.DataFrame:
    """Legacy function for polars output."""
    return fetch_wbp_biomart(xml_query, engine="polars", timeout=timeout)

# Convenience function for populateGeneList.py
def fetch_biomart_with_validation(
    xml_query: str, 
    expected_columns: set,
    timeout: int = REQUEST_TIMEOUT
) -> pl.DataFrame:
    """
    Enhanced fetch specifically for populateGeneList.py with column validation.
    """
    df = fetch_wbp_biomart(xml_query, engine="polars", timeout=timeout)
    
    # Immediate column validation (aligns with populateGeneList.py)
    _validate_biomart_columns(df, expected_columns)
    
    return df

def _validate_biomart_columns(df: pl.DataFrame, expected_columns: set) -> None:
    """Validate that BioMart returned expected columns."""
    if df.is_empty():
        raise ValueError("BioMart returned empty DataFrame")
    
    actual_columns = set(df.columns)
    missing_columns = expected_columns - actual_columns
    
    if missing_columns:
        available_columns = ", ".join(sorted(actual_columns))
        raise ValueError(
            f"BioMart response missing columns: {sorted(missing_columns)}\n"
            f"Available: {available_columns}"
        )


