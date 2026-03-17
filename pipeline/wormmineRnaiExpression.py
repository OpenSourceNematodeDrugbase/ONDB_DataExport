"""
wormmineRnaiExpression.py

Python 3.12–compatible WormMine expression fetcher.

Step 1: single-gene expression (long format)
Step 2: multi-gene expression with batched parallel fetching (strategy C)
Step 3: aggregation per C. elegans gene (one row per CE gene)

Single-gene (long):
    fetch_wormmine_expression_for_gene(gene_id)

Multi-gene (long, many rows per gene):
    fetch_wormmine_expression_for_genes(gene_ids, batch_size=20, max_workers=4)

Aggregated (one row per CE gene):
    aggregate_expression_by_gene(df_long)
"""

from __future__ import annotations

from typing import Iterable, List, Optional
import concurrent.futures

import requests
import pandas as pd

# WormMine endpoint (same as RNAi module; http to match UI behaviour)
WORMMINE_URL = "http://intermine.wormbase.org/tools/wormmine/service/query/results"
DEFAULT_TIMEOUT = 60  # seconds

# Long-format columns for expression-only query
EXPRESSION_COLUMNS = [
    # gene identifiers
    "wormbase_gene_id",              # Gene.primaryIdentifier
    "sequence_name",                 # Gene.secondaryIdentifier
    "gene_name",                     # Gene.symbol

    # expression patterns (gene-level)
    "expr_remark",                   # Gene.expressionPatterns.remark
    "expr_reporter_gene",            # Gene.expressionPatterns.reporterGene
    "expr_primary_id",               # Gene.expressionPatterns.primaryIdentifier
    "expr_pattern",                  # Gene.expressionPatterns.pattern
    "expr_subcellular_location",     # Gene.expressionPatterns.subcellularLocalization

    # expression patterns nested under lifeStages
    "expr_ls_remark",                # Gene.expressionPatterns.lifeStages.expressionPatterns.remark
    "expr_ls_reporter_gene",         # Gene.expressionPatterns.lifeStages.expressionPatterns.reporterGene
    "expr_ls_primary_id",            # Gene.expressionPatterns.lifeStages.expressionPatterns.primaryIdentifier
    "expr_ls_pattern",               # Gene.expressionPatterns.lifeStages.expressionPatterns.pattern
    "expr_ls_subcellular_location",  # Gene.expressionPatterns.lifeStages.expressionPatterns.subcellularLocalization
]


def _build_expression_query_xml(raw_gene_id: str) -> str:
    """
    Build the WormMine XML query for expression-only data for a single gene.

    This is derived from the extended RNAi+expression template, but restricted
    to the Gene + expression / life-stage expression paths only.
    """
    gene_id = raw_gene_id.strip()

    xml = f"""<query model="genomic"
  view="Gene.primaryIdentifier Gene.secondaryIdentifier Gene.symbol Gene.expressionPatterns.remark Gene.expressionPatterns.reporterGene Gene.expressionPatterns.primaryIdentifier Gene.expressionPatterns.pattern Gene.expressionPatterns.subcellularLocalization Gene.expressionPatterns.lifeStages.expressionPatterns.remark Gene.expressionPatterns.lifeStages.expressionPatterns.reporterGene Gene.expressionPatterns.lifeStages.expressionPatterns.primaryIdentifier Gene.expressionPatterns.lifeStages.expressionPatterns.pattern Gene.expressionPatterns.lifeStages.expressionPatterns.subcellularLocalization"
  sortOrder="Gene.primaryIdentifier ASC">
  <constraint path="Gene.primaryIdentifier" op="=" value="{gene_id}" code="A" />
</query>
"""
    return xml


def fetch_wormmine_expression_for_gene(gene_id: str) -> pd.DataFrame:
    """
    Fetch expression rows for a single WormBase gene (C. elegans).

    Returns a long-format DataFrame with one row per expression record
    (including life-stage expression), or an empty frame with the same
    columns if no data / error.
    """
    if not gene_id or not isinstance(gene_id, str):
        raise ValueError("gene_id must be a non-empty string")

    xml = _build_expression_query_xml(gene_id)

    try:
        response = requests.post(
            WORMMINE_URL,
            params={"format": "json"},
            data={"query": xml},
            timeout=DEFAULT_TIMEOUT,
        )
        response.raise_for_status()
    except Exception as exc:  # pragma: no cover - network failures
        print(f"[wormmine expression] HTTP error for {gene_id}: {exc}")
        return pd.DataFrame(columns=EXPRESSION_COLUMNS)

    try:
        data = response.json()
    except ValueError as exc:  # pragma: no cover - bad JSON
        print(f"[wormmine expression] JSON parse error for {gene_id}: {exc}")
        return pd.DataFrame(columns=EXPRESSION_COLUMNS)

    raw_rows = data.get("results", [])
    if not raw_rows:
        return pd.DataFrame(columns=EXPRESSION_COLUMNS)

    mapped: List[dict] = []

    for row in raw_rows:
        # Guard against unexpected row length
        if len(row) != len(EXPRESSION_COLUMNS):
            # Soft-fail: skip malformed rows but keep going
            print(
                f"[wormmine expression] Skipping row for {gene_id}: "
                f"expected {len(EXPRESSION_COLUMNS)} cols, got {len(row)}"
            )
            continue

        mapped.append(
            {
                "wormbase_gene_id": row[0],
                "sequence_name": row[1],
                "gene_name": row[2],
                "expr_remark": row[3],
                "expr_reporter_gene": row[4],
                "expr_primary_id": row[5],
                "expr_pattern": row[6],
                "expr_subcellular_location": row[7],
                "expr_ls_remark": row[8],
                "expr_ls_reporter_gene": row[9],
                "expr_ls_primary_id": row[10],
                "expr_ls_pattern": row[11],
                "expr_ls_subcellular_location": row[12],
            }
        )

    return pd.DataFrame(mapped, columns=EXPRESSION_COLUMNS)


# --- helper: chunking + strategy C multi-gene fetch ---


def _chunked(iterable: Iterable[str], size: int) -> Iterable[List[str]]:
    """Yield successive chunks of size `size` from `iterable`."""
    batch: List[str] = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def fetch_wormmine_expression_for_genes(
    gene_ids: Iterable[str],
    batch_size: int = 20,
    max_workers: int = 4,
) -> pd.DataFrame:
    """
    Fetch expression data for many C. elegans genes (long format).

    Strategy C:
      - clean & de-duplicate gene IDs
      - process in batches of `batch_size`
      - inside each batch, call single-gene fetch in parallel
      - concatenate all results (long format)
    """
    # Clean & deduplicate
    cleaned: List[str] = []
    seen: set[str] = set()
    for gid in gene_ids:
        if gid is None:
            continue
        s = str(gid).strip()
        if not s:
            continue
        if s in seen:
            continue
        seen.add(s)
        cleaned.append(s)

    if not cleaned:
        return pd.DataFrame(columns=EXPRESSION_COLUMNS)

    all_frames: List[pd.DataFrame] = []

    for batch in _chunked(cleaned, batch_size):
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(fetch_wormmine_expression_for_gene, gid): gid
                for gid in batch
            }

            for fut in concurrent.futures.as_completed(futures):
                gid = futures[fut]
                try:
                    df = fut.result()
                except Exception as exc:  # pragma: no cover - unexpected
                    print(f"[wormmine expression] Error fetching {gid}: {exc}")
                    continue

                if df is not None and not df.empty:
                    all_frames.append(df)

    if not all_frames:
        return pd.DataFrame(columns=EXPRESSION_COLUMNS)

    long_df = pd.concat(all_frames, ignore_index=True)
    return long_df


# --- aggregation to one row per gene ---


def _join_unique(values: Iterable[Optional[str]]) -> str:
    """Helper: collect unique, non-empty strings and join with '; '."""
    seen: list[str] = []
    for v in values:
        if v is None:
            continue
        s = str(v).strip()
        if not s:
            continue
        if s not in seen:
            seen.append(s)
    return "; ".join(seen)


def _compute_expr_flags_for_group(group: pd.DataFrame) -> dict:
    """
    Compute expression + life-stage flags for a single gene group.

    Uses the long-format columns:
      - expr_pattern
      - expr_subcellular_location
      - expr_ls_pattern
      - expr_ls_subcellular_location
    """
    text_bits: list[str] = []

    expr_cols = [
        "expr_pattern",
        "expr_subcellular_location",
        "expr_ls_pattern",
        "expr_ls_subcellular_location",
    ]
    for col in expr_cols:
        if col not in group.columns:
            continue
        vals = group[col].dropna().astype(str)
        if not vals.empty:
            text_bits.extend(vals.tolist())

    full_text = " ".join(t.lower() for t in text_bits if t)
    full_text = full_text.strip()

    if not full_text:
        # No expression information at all
        return {key: False for key in GENE_EXPR_FLAG_COLUMNS}

    def _contains_any(keywords: list[str]) -> bool:
        return any(kw in full_text for kw in keywords)

    # Simple presence flag
    any_expression = bool(full_text)

    # --- tissues / systems ---
    neuronal_kw = [
        "neuron", "neurons", "neuronal", "nerve ring", "nerve cord",
        "ganglion", "ganglia", "sensory neuron", "motor neuron",
    ]
    muscle_kw = ["muscle", "body wall muscle", "pharyngeal muscle", "vulval muscle"]
    intestinal_kw = ["intestine", "intestinal", "gut"]
    hypodermis_kw = ["hypodermis", "hypodermal", "skin", "epidermis"]
    pharynx_kw = ["pharynx", "pharyngeal"]
    germline_kw = [
        "germline", "germ line", "gonad", "gonadal", "gonad arm",
        "oocyte", "oocytes", "spermatheca", "reproductive system", "uterus",
    ]
    coelomocyte_kw = ["coelomocyte", "coelomocytes"]

    # --- subcellular keywords ---
    nuclear_kw = ["nuclear", "nucleus", "nuclei"]
    cytoplasmic_kw = ["cytoplasm", "cytoplasmic"]
    membrane_kw = ["membrane", "cell membrane", "plasma membrane", "cell surface"]

    # --- life-stage keywords (mainly from life-stage pattern text) ---
    embryo_kw = [
        "embryo", "embryonic", "precomma", "comma stage", "postcomma",
        "150-cell", "200-cell",
    ]
    larval_kw = [
        "larva", "larval", "larvae", " l1", " l2", " l3", " l4",
    ]
    adult_kw = [
        "adult", "young adult", "adult hermaphrodite", "adult male",
    ]
    male_specific_kw = [
        "male tail", "adult male", "male-specific", "male specific",
    ]

    flags = {
        # presence
        "ce_expr_any_expression": any_expression,

        # tissues / systems
        "ce_expr_any_neuronal": _contains_any(neuronal_kw),
        "ce_expr_any_muscle": _contains_any(muscle_kw),
        "ce_expr_any_intestinal": _contains_any(intestinal_kw),
        "ce_expr_any_hypodermis": _contains_any(hypodermis_kw),
        "ce_expr_any_pharynx": _contains_any(pharynx_kw),
        "ce_expr_any_germline": _contains_any(germline_kw),
        "ce_expr_any_coelomocyte": _contains_any(coelomocyte_kw),

        # subcellular
        "ce_expr_any_nuclear": _contains_any(nuclear_kw),
        "ce_expr_any_cytoplasmic": _contains_any(cytoplasmic_kw),
        "ce_expr_any_membrane": _contains_any(membrane_kw),

        # life stages
        "ce_expr_any_embryo": _contains_any(embryo_kw),
        "ce_expr_any_larval": _contains_any(larval_kw),
        "ce_expr_any_adult": _contains_any(adult_kw),
        "ce_expr_any_male_specific": _contains_any(male_specific_kw),
    }

    return flags


GENE_EXPR_FLAG_COLUMNS = [
    # presence
    "ce_expr_any_expression",

    # tissues / systems
    "ce_expr_any_neuronal",
    "ce_expr_any_muscle",
    "ce_expr_any_intestinal",
    "ce_expr_any_hypodermis",
    "ce_expr_any_pharynx",
    "ce_expr_any_germline",
    "ce_expr_any_coelomocyte",

    # subcellular
    "ce_expr_any_nuclear",
    "ce_expr_any_cytoplasmic",
    "ce_expr_any_membrane",

    # life stages
    "ce_expr_any_embryo",
    "ce_expr_any_larval",
    "ce_expr_any_adult",
    "ce_expr_any_male_specific",
]


GENE_ID_COLUMNS = ["wormbase_gene_id", "sequence_name", "gene_name"]

GENE_OUTPUT_COLUMNS_EXPR = GENE_ID_COLUMNS + GENE_EXPR_FLAG_COLUMNS


def aggregate_expression_by_gene(df_long: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate long-format expression rows to one row per gene.

    Parameters
    ----------
    df_long : pandas.DataFrame
        Output of fetch_wormmine_expression_for_genes (or _for_gene),
        one row per expression record.

    Returns
    -------
    pandas.DataFrame
        One row per wormbase_gene_id with expression flags.
    """
    if df_long is None or df_long.empty:
        return pd.DataFrame(columns=GENE_OUTPUT_COLUMNS_EXPR)

    required_cols = ["wormbase_gene_id", "sequence_name", "gene_name"]
    missing = [c for c in required_cols if c not in df_long.columns]
    if missing:
        raise ValueError(f"aggregate_expression_by_gene: missing columns {missing}")

    grouped = df_long.groupby("wormbase_gene_id", sort=True, dropna=False)

    records: list[dict] = []

    for gid, group in grouped:
        record: dict = {}

        record["wormbase_gene_id"] = gid
        record["sequence_name"] = _join_unique(group["sequence_name"])
        record["gene_name"] = _join_unique(group["gene_name"])

        # Expression flags
        expr_flags = _compute_expr_flags_for_group(group)
        record.update(expr_flags)

        records.append(record)

    agg_df = pd.DataFrame.from_records(records, columns=GENE_OUTPUT_COLUMNS_EXPR)

    # Assert one row per gene ID
    if agg_df["wormbase_gene_id"].duplicated().any():
        raise AssertionError(
            "aggregate_expression_by_gene: duplicate wormbase_gene_id in output"
        )

    return agg_df


def get_gene_expression_summary(
    gene_ids: Iterable[str],
    batch_size: int = 20,
    max_workers: int = 4,
) -> pd.DataFrame:
    """
    Convenience wrapper: fetch long-format expression for many genes
    and return the aggregated, one-row-per-gene summary.

    This is the function that ONDB's run_v01.py is likely to call.
    """
    df_long = fetch_wormmine_expression_for_genes(
        gene_ids,
        batch_size=batch_size,
        max_workers=max_workers,
    )
    return aggregate_expression_by_gene(df_long)
