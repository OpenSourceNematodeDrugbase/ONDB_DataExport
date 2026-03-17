# run_v01.py — Phase 1A: per-run folder + logging skeleton (no BioMart yet)


# make local helper modules importable (same approach used in runPipeline.py)
import json, sys, os
from datetime import datetime
from pathlib import Path
import pandas as pd
from populateGeneList_2 import *

# ---- CONFIG (adjust if needed) ----


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

OUTDIR = Path("pipeline_runs")   # was Path("/pipeline")
LABEL = "biomart-v01"
SPECIES_LIST = ["trtricprjeb535"]  # trial set

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


# # --- Phase 2: load canonical schema profile once per run ---


# SCHEMA_PATH = Path(__file__).parent  / "config" / "schema_profiles.json"


# try:
#     with SCHEMA_PATH.open("r", encoding="utf-8") as f:
#         SCHEMA_PROFILES = json.load(f)
# except FileNotFoundError:
#     raise RuntimeError(f"Schema profile not found at {SCHEMA_PATH}. Check exists?")

# COLUMNS_SPEC_VERSION = SCHEMA_PROFILES.get("columns_spec_version", "unknown")
# GENES_PROFILE = SCHEMA_PROFILES["profiles"]["GENES_BASE"]
# CDS_PROFILE   = SCHEMA_PROFILES["profiles"]["CDS_BASE"]

# print(f"[schema] Loaded {COLUMNS_SPEC_VERSION} from {SCHEMA_PATH}")




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
    Minimal per-species fetch using your current helper(s).
    Goal for this step: return a DataFrame with at least the gene IDs.
    """
    print(f"[info] Fetching base gene list for species: {species_code}")
    # Your helper likely accepts a species code string (you used GENOMES before).
    # If your function name differs (e.g., getGeneList), just swap the call below.
    gene_df = retrieveGeneListFromWbpBiomart(species_code)

    # Ensure it's a DataFrame
    if not isinstance(gene_df, pd.DataFrame):
        gene_df = pd.DataFrame(gene_df)

    # Add a clear species column so we can concatenate cleanly
    gene_df["species"] = species_code
    return gene_df





def main():
    run_dir, log_path = start_run_folder()
    sys.stdout = Tee(log_path)   # everything you print is also saved
    print(f"[info] Run folder: {run_dir}")
    print(f"[info] Species list: {SPECIES_LIST}")

    # # placeholder dataframe so the script produces real files in 1A
    df = pd.DataFrame({"species": SPECIES_LIST})
    df.to_csv(run_dir / "export.csv", index=False)
    (run_dir / "export.json").write_text(df.to_json(orient="records", indent=2), encoding="utf-8")
    write_run_json(run_dir, SPECIES_LIST)


        # --- replace the placeholder with this block ---
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

  

    meta_path = run_dir / "run.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))

    meta["counts"] = {
        "total_rows": int(combined.shape[0]),
        "by_species": {sp: int((combined["species"] == sp).sum()) for sp in SPECIES_LIST}
    }
   

    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")


    # Save outputs
    combined.to_csv(run_dir / "export.csv", index=False)
    (run_dir / "export.json").write_text(combined.to_json(orient="records"), encoding="utf-8")

   


if __name__ == "__main__":
    main()
