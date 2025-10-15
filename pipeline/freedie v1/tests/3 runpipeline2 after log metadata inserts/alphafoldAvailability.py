# alphafoldAvailability.py
import os, requests, pandas as pd
from io import StringIO

import os, sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)


# Optional: only needed if you want global metadata logging (Option B)
try:
    from utils_logging import log_append, file_sha256
except Exception:
    log_append = None
    file_sha256 = None

MARTSERVICE = "https://parasite.wormbase.org/biomart/martservice"

def _build_xml(species_keys: str, excluded_flag: int) -> str:
    """
    excluded_flag = 0 -> genes WITH AlphaFold
    excluded_flag = 1 -> genes WITHOUT AlphaFold
    species_keys: single key or comma-separated list (e.g., "wubancprjeb536,wubancprjna275548")
    """
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="parasite_mart" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
  <Dataset name="wbps_gene" interface="default">
    <Filter name="with_alphafold" excluded="{excluded_flag}"/>
    <Filter name="species_id_1010" value="{species_keys}"/>
    <Attribute name="production_name_1010" />
    <Attribute name="wbps_gene_id" />
  </Dataset>
</Query>
"""

def _post(xml_query: str):
    r = requests.post(MARTSERVICE, data={"query": xml_query}, timeout=60)
    r.raise_for_status()
    return r.text, {k.lower(): v for k, v in r.headers.items()}

def _query_block(species_keys: str, excluded_flag: int) -> tuple[pd.DataFrame, str, str, dict]:
    """
    Returns:
      df (production_name, Gene stable ID), xml_string, tsv_string, response_headers
    """
    xml = _build_xml(species_keys, excluded_flag)
    tsv, headers = _post(xml)
    if not tsv.strip():
        df = pd.DataFrame(columns=["production_name", "Gene stable ID"])
    else:
        df = pd.read_csv(StringIO(tsv), sep="\t", header=None, names=["production_name", "Gene stable ID"])
        df = df.drop_duplicates()
    return df, xml, tsv, headers

def queryAlphaFoldAvailability(
    genomes: str,
    log: bool = True,
    logs_dir: str = "logs/alphafold",
    meta_path: str = "logs/RUN_METADATA.json"
) -> pd.DataFrame:
    """
    Input:
        genomes: WBPS species key(s), comma-separated.
    Output:
        DataFrame with columns: Gene stable ID, production_name, alphafold_available
    Side effect (if log=True and utils_logging is available):
        - Saves XML/TSV for WITH and WITHOUT queries
        - Appends an 'alphafold' stage to pipeline/logs/RUN_METADATA.json
    """
    # Run the two queries
    df_with, xml_with, tsv_with, hdr_with = _query_block(genomes, excluded_flag=0)
    df_without, xml_without, tsv_without, hdr_without = _query_block(genomes, excluded_flag=1)

    # Merge with True winning on duplicates (defensive)
    df_with["alphafold_available"] = True
    df_without["alphafold_available"] = False
    all_df = pd.concat([df_with, df_without], ignore_index=True)
    all_df = (all_df.sort_values("alphafold_available", ascending=False)
                    .drop_duplicates(subset=["Gene stable ID"], keep="first"))
    all_df = all_df[["Gene stable ID", "production_name", "alphafold_available"]].sort_values("Gene stable ID")

    # Quick console stats
    genes_total = all_df["Gene stable ID"].nunique()
    genes_with = int((all_df["alphafold_available"] == True).sum())
    genes_without = genes_total - genes_with
    genomes_total = all_df["production_name"].nunique()
    print(f"[AlphaFold] GENES -> With: {genes_with:,} | Without: {genes_without:,} | Total: {genes_total:,}")
    print(f"[AlphaFold] GENOMES (distinct production_name): {genomes_total:,}")

    # Optional logging: write XML/TSV and append to global metadata
    if log and log_append is not None:
        os.makedirs(logs_dir, exist_ok=True)

        xml_with_path = os.path.join(logs_dir, "alphafold_with.xml")
        xml_without_path = os.path.join(logs_dir, "alphafold_without.xml")
        tsv_with_path = os.path.join(logs_dir, "alphafold_with.tsv")
        tsv_without_path = os.path.join(logs_dir, "alphafold_without.tsv")
        out_csv_path = "alphafold_availability.csv"

        # Save raw artifacts
        with open(xml_with_path, "w", encoding="utf-8") as f: f.write(xml_with)
        with open(xml_without_path, "w", encoding="utf-8") as f: f.write(xml_without)
        with open(tsv_with_path, "w", encoding="utf-8") as f: f.write(tsv_with)
        with open(tsv_without_path, "w", encoding="utf-8") as f: f.write(tsv_without)

        # Append a stage record into the global metadata file
        log_append("alphafold", {
            "genomes": genomes,
            "queries": {
                "with_alphafold": {
                    "xml_file": xml_with_path,
                    "tsv_file": tsv_with_path,
                    "rows": int(df_with.shape[0]),
                    "xml_sha256": file_sha256(xml_with_path) if file_sha256 else None,
                    "tsv_sha256": file_sha256(tsv_with_path) if file_sha256 else None
                },
                "without_alphafold": {
                    "xml_file": xml_without_path,
                    "tsv_file": tsv_without_path,
                    "rows": int(df_without.shape[0]),
                    "xml_sha256": file_sha256(xml_without_path) if file_sha256 else None,
                    "tsv_sha256": file_sha256(tsv_without_path) if file_sha256 else None
                }
            },
            "counts": {
                "genes_total": genes_total,
                "genes_with": genes_with,
                "genes_without": genes_without,
                "genomes_total": genomes_total
            },
            "outputs": {
                "csv_file": out_csv_path
            }
        }, meta_path=meta_path)
    

    # # Optional logging: always create folder & write artifacts if log=True.
    # # If utils_logging is available, also append to the global RUN_METADATA.json.
    # if log:
    #     os.makedirs(logs_dir, exist_ok=True)

    #     xml_with_path = os.path.join(logs_dir, "alphafold_with.xml")
    #     xml_without_path = os.path.join(logs_dir, "alphafold_without.xml")
    #     tsv_with_path = os.path.join(logs_dir, "alphafold_with.tsv")
    #     tsv_without_path = os.path.join(logs_dir, "alphafold_without.tsv")
    #     out_csv_path = "alphafold_availability.csv"

    #     # Save raw artifacts regardless of utils_logging presence
    #     with open(xml_with_path, "w", encoding="utf-8") as f: f.write(xml_with)
    #     with open(xml_without_path, "w", encoding="utf-8") as f: f.write(xml_without)
    #     with open(tsv_with_path, "w", encoding="utf-8") as f: f.write(tsv_with)
    #     with open(tsv_without_path, "w", encoding="utf-8") as f: f.write(tsv_without)

    #     # Append a stage record only if utils_logging is available
    #     if log_append is not None:
    #         log_append("alphafold", {
    #             "genomes": genomes,
    #             "queries": {
    #                 "with_alphafold": {
    #                     "xml_file": xml_with_path,
    #                     "tsv_file": tsv_with_path,
    #                     "rows": int(df_with.shape[0]),
    #                     "xml_sha256": file_sha256(xml_with_path) if file_sha256 else None,
    #                     "tsv_sha256": file_sha256(tsv_with_path) if file_sha256 else None
    #                 },
    #                 "without_alphafold": {
    #                     "xml_file": xml_without_path,
    #                     "tsv_file": tsv_without_path,
    #                     "rows": int(df_without.shape[0]),
    #                     "xml_sha256": file_sha256(xml_without_path) if file_sha256 else None,
    #                     "tsv_sha256": file_sha256(tsv_without_path) if file_sha256 else None
    #                 }
    #             },
    #             "counts": {
    #                 "genes_total": genes_total,
    #                 "genes_with": genes_with,
    #                 "genes_without": genes_without,
    #                 "genomes_total": genomes_total
    #             },
    #             "outputs": {
    #                 "csv_file": out_csv_path
    #             }
    #         }, meta_path=meta_path)




    return all_df
