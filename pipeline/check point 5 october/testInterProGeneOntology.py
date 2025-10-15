# testInterProGeneOntology.py — Step 2 (GO logic added)
# - EXACT same XML as before (header="1"); no further XML edits here
# - species_code single-species validation
# - Header discovery prints (already working)
# - Use interpro2go mapping + GO DAG (if available) to check target GO + descendants
# - Output: ["Gene stable ID", <test_name> (bool), <test_name>_evidence (str)]
# - Polars internally; Pandas at boundary
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Dict, Set, Iterable

import polars as pl
import pandas as pd

# Central fetcher (use whichever name exists in your repo)
try:
    from queryWbpBiomart import fetch_wbp_biomart_using_xml_polars
except ImportError:
    from queeryingwbpbiomart import fetch_wbp_biomart_using_xml_polars





def _validate_single_species(species_code: str) -> str:
    if not isinstance(species_code, str) or not species_code.strip():
        raise ValueError("[error] species_code must be a non-empty string")
    if "," in species_code:
        raise ValueError(f"[error] Exactly one species_code expected (no commas): {species_code!r}")
    return species_code.strip()


def _build_xml(species_code: str) -> str:
    # EXACT XML provided by you. No edits beyond species_code interpolation.
    return f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
              
        <Dataset name = "wbps_gene" interface = "default" >
        <Filter name = "species_id_1010" value = "{species_code}"/>
        <Filter name = "biotype" value = "protein_coding"/>
        <Attribute name = "wbps_gene_id" />
        <Attribute name = "interpro_id" />
        <Attribute name = "interpro_description" />
        </Dataset>
    </Query>
    """


# ---------- Helpers: interpro2go + GO descendants (no network calls) ----------

def _local_path(filename: str) -> Path:
    """Return a path to a file that sits next to this module."""
    return Path(__file__).parent / filename

def _load_interpro2go_map() -> Dict[str, Set[str]]:
    """
    Load 'interpro2go' from the SAME DIRECTORY as this script.
    Returns: dict like {"IPR000001": {"GO:0000001","GO:0000002"}, ...}
    If the file is missing, returns an empty dict (module will still run).
    """
    path = _local_path("interpro2go")
    ipr2go: Dict[str, Set[str]] = {}
    if not path.exists():
        print("[warn] interpro2go not found in this folder; proceeding with empty mapping.")
        return ipr2go

    import re
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            m_ipr = re.search(r"InterPro:(IPR\d+)", line)
            if not m_ipr:
                continue
            ipr = m_ipr.group(1)
            gos = set(re.findall(r"(GO:\d+)", line))
            if not gos:
                continue
            ipr2go.setdefault(ipr, set()).update(gos)

    print(f"[info] Loaded interpro2go from {path} "
          f"({sum(len(v) for v in ipr2go.values())} IPR→GO links across {len(ipr2go)} IPRs).")
    return ipr2go

def _load_go_children() -> Dict[str, Set[str]]:
    """
    Build parent->children map using ONLY 'is_a:' edges from 'go-basic.obo'
    located in the SAME DIRECTORY as this script.
    """
    path = _local_path("go-basic.obo")
    parent_to_children: Dict[str, Set[str]] = {}

    if not path.exists():
        print("[warn] go-basic.obo not found in this folder; GO descendants will NOT be expanded.")
        return parent_to_children

    current_id: str | None = None
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for raw in fh:
                line = raw.strip()
                if line == "[Term]":
                    current_id = None
                    continue
                if line.startswith("id: GO:"):
                    current_id = line.split("id: ")[1].strip()
                    continue
                if current_id is None:
                    continue
                if line.startswith("is_a: GO:"):
                    parent = line.split("is_a: ")[1].split()[0].strip()
                    parent_to_children.setdefault(parent, set()).add(current_id)
        print(f"[info] Loaded GO DAG (is_a only) from {path} "
              f"with {len(parent_to_children)} parent nodes.")
    except Exception as e:
        print(f"[warn] Failed reading {path}: {e}. Descendants will not expand.")
        return {}

    return parent_to_children

def _descendants(seed: str, parent_to_children: Dict[str, Set[str]]) -> Set[str]:
    """All descendants of seed including multi-level, not including the seed itself."""
    seen: Set[str] = set()
    stack = [seed]
    while stack:
        node = stack.pop()
        for child in parent_to_children.get(node, ()):
            if child not in seen:
                seen.add(child)
                stack.append(child)
    return seen


# ---------- Main ----------
def testInterProGeneOntology(species_code: str, go_term: str, test_name: str) -> pd.DataFrame:
    """
    Returns a Pandas DataFrame with columns:
      - "Gene stable ID"
      - <test_name> (bool)
      - <test_name>_evidence (str)

    Notes:
      * Requires Biomart to return columns: "Gene stable ID", "InterPro ID", "InterPro description".
      * If interpro2go/go-basic.obo are missing, descendants may not expand and mapping may be empty.
      * Genes with no InterPro rows are not emitted here; they appear after your final left-join.
    """
    species_code = _validate_single_species(species_code)
    xml_query = _build_xml(species_code)

    # Fetch
    df: pl.DataFrame = fetch_wbp_biomart_using_xml_polars(xml_query)

    # Header discovery
    cols = df.columns
    print(f"[info] Biomart header discovered ({len(cols)} columns):")
    for i, c in enumerate(cols, start=1):
        print(f"  [{i:02d}] {c}")

    # Required columns
    required = {"Gene stable ID", "InterPro ID", "InterPro description"}
    missing = [c for c in required if c not in cols]
    if missing:
        raise RuntimeError(f"[error] Missing required columns: {missing}. Check XML and Biomart headers.")

    # Load mappings / DAG (best-effort; non-fatal if absent)
    ipr2go = _load_interpro2go_map()
    go_children = _load_go_children()

    # Allowed GO set = target + descendants (if DAG present)
    allowed_go: Set[str] = {go_term}
    if go_children:
        allowed_go |= _descendants(go_term, go_children)

    # Prepare row-level flags and labels
    df2 = (
        df
        .filter(pl.col("InterPro ID").is_not_null() & (pl.col("InterPro ID") != ""))
        .with_columns([
            pl.struct(["InterPro ID", "InterPro description"]).map_elements(
                lambda s: f"{s['InterPro ID']} - {s['InterPro description']}"
            ).alias("_ipr_label"),
            pl.col("InterPro ID").map_elements(
                lambda ipr: bool(ipr2go.get(ipr, set()) & allowed_go)
            ).alias("_ipr_matches"),
        ])
    )

    if df2.height == 0:
        # No InterPro rows at all for this species; return empty shell with right schema
        empty = pl.DataFrame(
            {
                "Gene stable ID": pl.Series([], dtype=pl.Utf8),
                test_name: pl.Series([], dtype=pl.Boolean),
                f"{test_name}_evidence": pl.Series([], dtype=pl.Utf8),
            }
        )
        return empty.to_pandas()

    # Aggregate per gene
    grouped = (
        df2
        .group_by("Gene stable ID", maintain_order=True)
        .agg([
            pl.any("_ipr_matches").alias(test_name),
            # collect only matching labels into a list
            pl.col("_ipr_label").filter(pl.col("_ipr_matches")).alias("_match_labels"),
        ])
        .with_columns([
            pl.when(pl.col(test_name))
              .then(
                  pl.lit(f"Encodes protein with InterPro domain(s) that are related to {go_term}: ")
                  + pl.col("_match_labels").list.unique().list.sort().list.join("; ")
              )
              .otherwise(pl.lit(f"Does not have InterPro domains related to {go_term}"))
              .alias(f"{test_name}_evidence"),
        ])
        .select(["Gene stable ID", test_name, f"{test_name}_evidence"])
    )

    # Integrity: one row per gene
    n = grouped.height
    n_unique = grouped.select(pl.col("Gene stable ID")).unique().height
    if n != n_unique:
        raise AssertionError("[error] Duplicate 'Gene stable ID' rows detected in module output.")

    return grouped.to_pandas()
