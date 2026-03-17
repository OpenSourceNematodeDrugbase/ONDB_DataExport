"""
wormmineRnaiPhenotypes.py

Python 3.12–compatible WormMine RNAi phenotype fetcher.

Step 1: single-gene RNAi (long format)
Step 2: multi-gene RNAi with batched parallel fetching (strategy C)
Step 3: aggregation per C. elegans gene (one row per CE gene)
----------------------------------------------------------------

Single-gene (long):
    fetch_wormmine_rnai_for_gene(gene_id)

Multi-gene (long, many rows per gene):
    fetch_wormmine_rnai_for_genes(gene_ids, batch_size=20, max_workers=4)

Aggregated (one row per CE gene):
    aggregate_rnai_by_gene(df_long)
"""

from __future__ import annotations

from typing import List, Any, Iterable, Optional
import concurrent.futures

import requests
import pandas as pd


# WormMine endpoint (REST API, http to match UI behaviour)
WORMMINE_URL = "http://intermine.wormbase.org/tools/wormmine/service/query/results"
DEFAULT_TIMEOUT = 60  # seconds

RNAI_COLUMNS = [
    "wormbase_gene_id",
    "sequence_name",
    "rnai_result_id",
    "phenotype_id",
    "phenotype_name",
    "phenotype_remark",
    "rnai_result_remark",
]


# ---------------------------------------------------------
#  Build the query XML (exactly as in the working template)
# ---------------------------------------------------------

def _build_rnai_query_xml(raw_gene_id: str) -> str:
    """
    Build the WormMine XML query dynamically for a single gene.

    Based on the working XML from WormMine QueryBuilder:
    <query model="genomic"
           view="Gene.primaryIdentifier Gene.secondaryIdentifier ..."
           sortOrder="Gene.primaryIdentifier ASC">
      <constraint path="Gene.primaryIdentifier" op="=" value="WBGene00006757" code="A" />
    </query>
    """
    gene_id = raw_gene_id.strip()  # avoid leading/trailing spaces

    xml = f"""<query model="genomic"
  view="Gene.primaryIdentifier Gene.secondaryIdentifier Gene.RNAiResult.primaryIdentifier Gene.RNAiResult.phenotype.identifier Gene.RNAiResult.phenotype.name Gene.RNAiResult.phenotypeRemark Gene.RNAiResult.remark"
  sortOrder="Gene.primaryIdentifier ASC">
  <constraint path="Gene.primaryIdentifier" op="=" value="{gene_id}" code="A" />
</query>
"""
    return xml


# ---------------------------------------------------------
#  Core: single-gene fetch (already tested)
# ---------------------------------------------------------

def fetch_wormmine_rnai_for_gene(gene_id: str) -> pd.DataFrame:
    """
    Fetch RNAi phenotype rows for a single WormBase gene.

    Parameters
    ----------
    gene_id : str
        WormBase Gene ID, e.g. "WBGene00006757".

    Returns
    -------
    pandas.DataFrame
        One row per (gene, RNAi result, phenotype), with columns:
            - wormbase_gene_id
            - sequence_name
            - rnai_result_id
            - phenotype_id
            - phenotype_name
            - phenotype_remark
            - rnai_result_remark

        If the gene has no RNAi phenotypes OR the query returns no rows,
        an empty DataFrame with the same columns is returned.
    """
    if not gene_id or not isinstance(gene_id, str):
        raise ValueError("gene_id must be a non-empty string")

    xml = _build_rnai_query_xml(gene_id)

    response = requests.post(
        WORMMINE_URL,
        params={"format": "json"},
        data={"query": xml},
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        timeout=DEFAULT_TIMEOUT,
    )

    try:
        response.raise_for_status()
    except requests.HTTPError as e:
        print(f"[error] WormMine HTTP error for {gene_id}:", e)
        print("[error] Response text (first 300 chars):")
        print(response.text[:300])
        # Return empty shell rather than killing the whole pipeline
        return pd.DataFrame(columns=RNAI_COLUMNS)

    try:
        data: Any = response.json()
    except ValueError:
        print(f"[error] WormMine response was not valid JSON for {gene_id}:")
        print(response.text[:300])
        return pd.DataFrame(columns=RNAI_COLUMNS)

    results = data.get("results", [])
    if not isinstance(results, list) or len(results) == 0:
        # No rows → return empty schema
        return pd.DataFrame(columns=RNAI_COLUMNS)

    mapped: List[dict] = []
    for row in results:
        if not isinstance(row, list) or len(row) < 7:
            continue  # skip malformed rows

        mapped.append(
            {
                "wormbase_gene_id": row[0],   # Gene.primaryIdentifier
                "sequence_name": row[1],      # Gene.secondaryIdentifier
                "rnai_result_id": row[2],     # RNAiResult.primaryIdentifier
                "phenotype_id": row[3],       # phenotype.identifier
                "phenotype_name": row[4],     # phenotype.name
                "phenotype_remark": row[5],   # phenotypeRemark
                "rnai_result_remark": row[6], # remark
            }
        )

    return pd.DataFrame(mapped, columns=RNAI_COLUMNS)


# ---------------------------------------------------------
#  Helper: chunking for batched parallel strategy (C)
# ---------------------------------------------------------

def _chunked(iterable: Iterable[str], size: int) -> Iterable[List[str]]:
    """
    Yield lists of up to `size` items from `iterable`.
    """
    batch: List[str] = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


# ---------------------------------------------------------
#  Multi-gene fetch (long format, batched + parallel)
# ---------------------------------------------------------

def fetch_wormmine_rnai_for_genes(
    gene_ids: Iterable[str],
    batch_size: int = 20,
    max_workers: int = 4,
) -> pd.DataFrame:
    """
    Fetch RNAi phenotypes for many C. elegans genes.

    Strategy C:
      - clean & de-duplicate gene IDs
      - process in batches of `batch_size`
      - inside each batch, call single-gene fetch in parallel
      - concatenate all results (long format)

    Parameters
    ----------
    gene_ids : Iterable[str]
        Collection of WormBase Gene IDs (e.g. from the CE orthologue table).
    batch_size : int, optional
        Number of genes per parallel batch (default 20).
    max_workers : int, optional
        Maximum number of worker threads per batch (default 4).

    Returns
    -------
    pandas.DataFrame
        Long-format table with the same columns as fetch_wormmine_rnai_for_gene().
        Multiple rows per gene if multiple RNAi phenotypes / experiments.
    """
    # Normalise & deduplicate gene IDs
    cleaned: List[str] = []
    seen = set()
    for gid in gene_ids:
        if not gid:
            continue
        if not isinstance(gid, str):
            gid = str(gid)
        gid = gid.strip()
        if not gid or gid in seen:
            continue
        seen.add(gid)
        cleaned.append(gid)

    if not cleaned:
        return pd.DataFrame(columns=RNAI_COLUMNS)

    all_frames: List[pd.DataFrame] = []

    for batch in _chunked(cleaned, batch_size):
        print(f"[info] WormMine RNAi — processing batch of {len(batch)} genes")

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_gid = {
                executor.submit(fetch_wormmine_rnai_for_gene, gid): gid
                for gid in batch
            }

            for future in concurrent.futures.as_completed(future_to_gid):
                gid = future_to_gid[future]
                try:
                    df = future.result()
                except Exception as e:
                    print(f"[warn] WormMine RNAi fetch failed for {gid}: {e}")
                    continue

                if not df.empty:
                    all_frames.append(df)

    if not all_frames:
        return pd.DataFrame(columns=RNAI_COLUMNS)

    combined = pd.concat(all_frames, ignore_index=True)

    # Optional: quick sanity prints
    print(
        f"[info] WormMine RNAi — combined long-format table: "
        f"{combined.shape[0]} rows for "
        f"{combined['wormbase_gene_id'].nunique()} genes"
    )
    return combined


# ---------------------------------------------------------
#  Aggregation: one row per CE gene
# ---------------------------------------------------------

def _join_unique(values: Iterable[Optional[str]]) -> str:
    """Helper: collect unique, non-empty strings and join with '; '."""
    seen = []
    for v in values:
        if v is None:
            continue
        s = str(v).strip()
        if not s:
            continue
        if s not in seen:
            seen.append(s)
    return "; ".join(seen)


def aggregate_rnai_by_gene(df_long: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate long-format RNAi table to one row per C. elegans gene.

    Parameters
    ----------
    df_long : pandas.DataFrame
        Output from fetch_wormmine_rnai_for_genes (long format).

    Returns
    -------
    pandas.DataFrame
        One row per wormbase_gene_id with columns:
          - wormbase_gene_id
          - sequence_name
          - ce_rnai_result_ids
          - ce_rnai_phenotype_ids
          - ce_rnai_phenotype_names
          - ce_rnai_phenotype_remarks
          - ce_rnai_rnai_result_remarks
          - ce_rnai_any_lethal      (bool)
          - ce_rnai_any_sterile     (bool)
          - ce_rnai_any_locomotion  (bool)
    """
    if df_long is None or df_long.empty:
        cols = [
            "wormbase_gene_id",
            "sequence_name",
            "ce_rnai_result_ids",
            "ce_rnai_phenotype_ids",
            "ce_rnai_phenotype_names",
            "ce_rnai_phenotype_remarks",
            "ce_rnai_rnai_result_remarks",
            "ce_rnai_any_lethal",
            "ce_rnai_any_sterile",
            "ce_rnai_any_locomotion",
        ]
        return pd.DataFrame(columns=cols)

    grouped = df_long.groupby("wormbase_gene_id", dropna=True, sort=True)

    records: List[dict] = []
    for gene_id, group in grouped:
        # sequence_name: first non-null, if any
        seq_name: Optional[str]
        non_null_seq = group["sequence_name"].dropna()
        seq_name = non_null_seq.iloc[0] if not non_null_seq.empty else None

        # Collect unique lists
        joined_result_ids = _join_unique(group["rnai_result_id"])
        joined_pheno_ids = _join_unique(group["phenotype_id"])
        joined_pheno_names = _join_unique(group["phenotype_name"])
        joined_pheno_remarks = _join_unique(group["phenotype_remark"])
        joined_rnai_remarks = _join_unique(group["rnai_result_remark"])

        # For flags, work on lowercase phenotype names + remarks
        text_blobs: List[str] = []
        for col in ("phenotype_name", "phenotype_remark", "rnai_result_remark"):
            for v in group[col]:
                if v is None:
                    continue
                s = str(v).strip()
                if s:
                    text_blobs.append(s.lower())

        any_lethal = any("lethal" in t for t in text_blobs)
        any_sterile = any("sterile" in t for t in text_blobs)
        any_locomotion = any(
            ("locomotion" in t) or ("movement" in t) for t in text_blobs
        )

        records.append(
            {
                "wormbase_gene_id": gene_id,
                "sequence_name": seq_name,
                "ce_rnai_result_ids": joined_result_ids,
                "ce_rnai_phenotype_ids": joined_pheno_ids,
                "ce_rnai_phenotype_names": joined_pheno_names,
                "ce_rnai_phenotype_remarks": joined_pheno_remarks,
                "ce_rnai_rnai_result_remarks": joined_rnai_remarks,
                "ce_rnai_any_lethal": bool(any_lethal),
                "ce_rnai_any_sterile": bool(any_sterile),
                "ce_rnai_any_locomotion": bool(any_locomotion),
            }
        )

    agg_df = pd.DataFrame.from_records(records)

    # Ensure boolean dtypes
    for col in ["ce_rnai_any_lethal", "ce_rnai_any_sterile", "ce_rnai_any_locomotion"]:
        agg_df[col] = agg_df[col].fillna(False).astype(bool)

    # Sort by gene ID for determinism
    agg_df = agg_df.sort_values("wormbase_gene_id").reset_index(drop=True)
    return agg_df


# ---------------------------------------------------------
#  Manual test
# ---------------------------------------------------------

if __name__ == "__main__":
    # Single-gene test (already working)
    test_gene = "WBGene00006757"
    df_one = fetch_wormmine_rnai_for_gene(test_gene)
    print("[single] head:")
    print(df_one.head())
    print(f"[single] Rows: {len(df_one)}")

    # Multi-gene test with 2–3 known genes
    test_genes = [
        "WBGene00006757",  # unc-26
        "WBGene00004803",  # dpy-11
        "WBGene00000912",  # rde-1
    ]
    df_multi = fetch_wormmine_rnai_for_genes(test_genes, batch_size=2, max_workers=2)
    print("[multi] head:")
    print(df_multi.head())
    print(
        f"[multi] Rows: {len(df_multi)}, "
        f"unique genes: {df_multi['wormbase_gene_id'].nunique() if not df_multi.empty else 0}"
    )

    # Aggregation test
    df_agg = aggregate_rnai_by_gene(df_multi)
    print("[agg] head:")
    print(df_agg.head())
    print(
        f"[agg] Rows (one per CE gene): {len(df_agg)}, "
        f"unique genes: {df_agg['wormbase_gene_id'].nunique() if not df_agg.empty else 0}"
    )
