# run_v01.py —


# make local helper modules importable (same approach used in runPipeline.py)
import json, sys, os
from datetime import datetime
from pathlib import Path
import pandas as pd
from populateGeneList import *
from wbpHumanOrthologues import (
    retrieveHumanOrthologuesFromWbpBiomart,
    EXPECTED_HUMAN_ORTHO_COLUMNS,   
)
from testInterProGeneOntology import testInterProGeneOntology
from testInterProStringSearch import testInterProStringSearch



# ---- CONFIG (adjust if needed) ----


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

OUTDIR = Path("pipeline_runs")  
LABEL = "biomart-v01"
SPECIES_LIST = ["trtricprjeb535","wubancprjna275548"]  # trial set

# ---- tiny logger to file + console ----
def start_run_folder():
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S") + f"_{LABEL}"
    run_dir = OUTDIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "log.txt"
    return run_dir, log_path

       

class Tee:
    """Write prints to both console and log file."""
    def __init__(self, log_path):
        self.terminal = sys.stdout
        self.log = open(log_path, "a", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

    def close(self):
        try:
            self.flush()
        finally:
            self.log.close()

    # (Optional) enables: with Tee(path) as tee:
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False







def write_run_json(run_dir, species_list):
    meta = {
        "run_id": run_dir.name,
        "label": LABEL,
        "created": datetime.now().isoformat(timespec="seconds"),
        "species_list": species_list,
        "counts": {}
    }
    (run_dir / "run.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")






def fetch_base_biomart_for_species(species_code: str) -> pd.DataFrame:


    """
    Per-species fetch:
      1) base gene list (principal table)
      2) human orthologues (new module)
      3) left-merge on 'Gene stable ID'
    Returns Pandas, ready to concat and export.
    """
    print(f"[info] Fetching base gene list for species: {species_code}")
    principal_df = retrieveGeneListFromWbpBiomart(species_code)  # Pandas
    print(f"[info] Principal gene rows for {species_code}: {len(principal_df):,}")


    print(f"[info] Fetching human orthologues for species: {species_code}")
    human_df = retrieveHumanOrthologuesFromWbpBiomart(species_code)  # Pandas
    print(f"[info] Human-orthologues rows for {species_code}: {len(human_df):,}")


    # Ensure join key is string on both sides
    principal_df["Gene stable ID"] = principal_df["Gene stable ID"].astype(str)
    human_df["Gene stable ID"] = human_df["Gene stable ID"].astype(str)

    # Merge (keep BioMart labels verbatim; no renaming)
    merged = principal_df.merge(
        human_df[
            [
                "Gene stable ID",
                "Human gene stable ID",
                "Human gene name",
                "Homology type",
                "% identity",
                "Human % identity",
                "lacks_WBP_human_orthologue",
                "lacks_WBP_human_orthologue_evidence",
                "best_WBP_human_orthologue_lt_40pct_identity",
                "best_WBP_human_orthologue_lt_40pct_identity_evidence",
            ]
        ],
        on="Gene stable ID",
        how="left",
    )

        # Merge sanity: rows should match principal (left-merge)
    print(f"[info] Post-merge rows for {species_code}: {merged.shape[0]:,} "
          f"(matches principal? {merged.shape[0] == principal_df.shape[0]})")

    # Ensure all expected BioMart columns from the human module made it into the export
    missing_in_export = sorted(set(EXPECTED_HUMAN_ORTHO_COLUMNS) - set(merged.columns))
    print(f"[info] Export includes all expected human-orthologue columns? {len(missing_in_export) == 0}")
    if missing_in_export:
        print(f"[warn] Missing expected human-orthologue columns in export for {species_code}: {missing_in_export}")

    # Ensure 'species' exists for your by-species counts later in run.json
    merged["species"] = species_code


    # ---- Module annotations (left-join onto the principal table) ----
    # NOTE: Modules return Pandas with Biomart headers preserved.
    # GO-based:
    df_is_enzyme = testInterProGeneOntology(species_code, "GO:0003824", "is_enzyme")
    df_is_kinase = testInterProGeneOntology(species_code, "GO:0004672", "is_kinase")

    # String-search (use original regex first; we can widen later if needed)
    df_is_gpcr   = testInterProStringSearch(
        species_code,
        r"G-protein coupled receptor|GPCR.[^k]",
        "is_gpcr",
    )

    for mod_df in [df_is_enzyme, df_is_kinase, df_is_gpcr]:
        merged = merged.merge(mod_df, on="Gene stable ID", how="left")

    # Fill booleans after left-join (rows with no InterPro stay False/empty)
    for col in ["is_enzyme", "is_kinase", "is_gpcr"]:
        if col in merged.columns:
            merged[col] = merged[col].infer_objects(copy=False).astype(bool)

    # Quick per-species sanity counts
    print(f"[info] {species_code} is_enzyme TRUE:", int(merged.get("is_enzyme", pd.Series([], dtype=bool)).sum()))
    print(f"[info] {species_code} is_kinase TRUE:", int(merged.get("is_kinase", pd.Series([], dtype=bool)).sum()))
    print(f"[info] {species_code} is_gpcr   TRUE:", int(merged.get("is_gpcr",   pd.Series([], dtype=bool)).sum()))

    dups = merged["Gene stable ID"].duplicated(keep=False)
    if dups.any():
        sample = merged.loc[dups, "Gene stable ID"].head(5).tolist()
        raise AssertionError(f"[error] Duplicate 'Gene stable ID' inside species {species_code}, e.g. {sample}")


    return merged



def main():
    run_dir, log_path = start_run_folder()

    previous_stdout = sys.stdout
    tee = Tee(log_path)
    sys.stdout = tee
    try:
        print(f"[info] Run folder: {run_dir}")
        print(f"[info] Species list: {SPECIES_LIST}")

        # Ensure run.json exists now that the early placeholder block is removed
        write_run_json(run_dir, SPECIES_LIST)

        # ---- real per-species fetch ----
        per_species = []
        for sp in SPECIES_LIST:
            try:
                df_sp = fetch_base_biomart_for_species(sp)
                per_species.append(df_sp)

            except Exception as e:
                print(f"[warn] Failed species {sp}: {e}")

        if len(per_species) == 0:
            print("[error] No species succeeded; aborting.")
            return


        combined = pd.concat(per_species, ignore_index=True)
        # prints how many columns we should have and rows
        print(f"[info] Combined export shape: {combined.shape[0]:,} rows × {combined.shape[1]:,} cols")
        print(f"[info] Final export column count: {combined.shape[1]}")
        


        # ---- update run.json counts ----
        meta_path = run_dir / "run.json"
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        meta["counts"] = {
            "total_rows": int(combined.shape[0]),
            "by_species": {sp: int((combined["species"] == sp).sum()) for sp in SPECIES_LIST},
        }
        meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")

        # 1) One row per gene (no duplicates)
        dup_counts = combined["Gene stable ID"].value_counts()
        dups = dup_counts[dup_counts > 1]
        if not dups.empty:
            sample = dups.index.tolist()[:5]
            more = len(dups) - len(sample)
            raise AssertionError(
                f"[error] Duplicate 'Gene stable ID' rows detected, e.g. {sample}"
                + (f" (+{more} more)" if more > 0 else "")
            )

        # 2) Ensure booleans are really bool (if these columns exist)
        for col in ["is_enzyme", "is_kinase", "is_gpcr"]:
            if col in combined.columns:
                combined[col] = combined[col].fillna(False).astype(bool)
                assert combined[col].dtype == bool, f"[error] {col} must be bool; got {combined[col].dtype}"

        # 3) Evidence column exists alongside each boolean
        for col in ["is_enzyme", "is_kinase", "is_gpcr"]:
            if col in combined.columns:
                ev = f"{col}_evidence"
                assert ev in combined.columns, f"[error] Missing evidence column: {ev}"

        # ---- final exports ----
        combined.to_csv(run_dir / "export.csv", index=False)
        (run_dir / "export.json").write_text(combined.to_json(orient="records"), encoding="utf-8")

    finally:
        # Always restore stdout and close the log file
        sys.stdout = previous_stdout
        tee.close()



if __name__ == "__main__":
    main()
