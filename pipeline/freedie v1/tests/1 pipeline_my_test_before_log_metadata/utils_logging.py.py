# utils_logging.py
"""
Tiny helper to maintain a single global metadata file:
  pipeline/logs/RUN_METADATA.json

Usage pattern:
  from utils_logging import log_append
  log_append("alphafold", {"queries": {...}, "outputs": {...}})

This will:
  - create the file on first use with a fresh run_id (UTC timestamp),
  - append a new stage record with a timestamp,
  - keep everything in one JSON for the run.
"""

import os, json, datetime, hashlib

def _ensure_dir(p):
    os.makedirs(os.path.dirname(p), exist_ok=True)

def file_sha256(path: str) -> str | None:
    if not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()

def _new_run_doc(run_id: str | None = None):
    return {
        "run_id": run_id or datetime.datetime.utcnow().strftime("%Y-%m-%dT%H-%M-%SZ"),
        "created_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "stages": []
    }

def log_append(stage_name: str, payload: dict,
               meta_path: str = "pipeline/logs/RUN_METADATA.json",
               run_id: str | None = None):
    """
    Append a stage entry into the global RUN_METADATA.json.
    If the file doesn't exist, it will be created.

    Parameters
    ----------
    stage_name : str
      e.g., "gene_list", "alphafold", "orthologues", "interpro"
    payload : dict
      Any JSON-serializable dictionary. Typical keys:
        - inputs / queries / outputs
        - counts
        - notes
    meta_path : str
      Where to store the metadata file.
    run_id : str | None
      Optional fixed run id (else a timestamp is used on first creation).
    """
    _ensure_dir(meta_path)
    if os.path.exists(meta_path):
        with open(meta_path, "r", encoding="utf-8") as f:
            doc = json.load(f)
    else:
        doc = _new_run_doc(run_id=run_id)

    entry = {"name": stage_name,
             "timestamp_utc": datetime.datetime.utcnow().isoformat() + "Z"}
    entry.update(payload or {})

    doc["stages"].append(entry)

    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=2)
