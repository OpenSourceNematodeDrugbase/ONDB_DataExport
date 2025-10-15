# ONDB BioMart Pipeline — Contract Version (headered XML, single-species modules)

This repository implements a **deterministic, headered, single-species** BioMart pipeline that assembles a clean export of parasite genes with annotations from InterPro and human orthologues. It follows a strict contract shared by all modules:

- **Headered BioMart output** (`header="1"`), human‑readable labels **verbatim**.
- **Central fetch** via `queryWbpBiomart.py`, with retries/timeout.
- **Single species input** (`species_code`), validated early (no commas).
- **Polars internally; Pandas at the boundary** (return values).
- **No renaming of BioMart labels** (we keep column names exactly as returned).
- **Deterministic row integrity** (assert no duplicate `"Gene stable ID"` at module and/or pipeline boundary).
- **Minimal, explicit required-column checks** after fetch.
- **Provenance & outputs** handled by `run_v01.py` (run folder, log, run.json).

> **Zero silent XML changes policy**: All BioMart XML templates are **identical to the originals** apart from two allowed adjustments: set `header="1"` and interpolate `{species_code}`. No other edits (attributes, order, whitespace) unless explicitly approved.


---

## 1) Repository layout (key files)

```
pipeline/
  populateGeneList_2.py         # principal gene list (per species; headered XML)
  queryWbpBiomart.py            # central BioMart client (retries, timeout; Polars/Pandas)
  testInterProGeneOntology.py   # InterPro→GO term test (GO + interpro2go; headered XML)
  testInterProStringSearch.py   # InterPro string/regex test (headered XML)
  run_v01.py                    # orchestrator: run folder, logging, joins, export
  interpro2go                   # REQUIRED mapping file (place here, same folder)
  go-basic.obo                  # REQUIRED GO DAG (place here, same folder)
```

> If you keep additional modules, put them next to these files and follow the same contract.


---

## 2) Python & dependencies

- Python **3.10+** (3.12 confirmed).
- Libraries:
  - `polars`
  - `pandas`
  - `requests`

Install:
```bash
pip install pandas polars requests
```

*(No `goatools` needed — the GO module is self-contained with a minimal OBO parse.)*


---

## 3) Required local data files

Place these **in the same folder** as the Python modules (e.g., `pipeline/`).

- **`interpro2go`** — InterPro ↔ GO mapping (plain text).  
  *Used by:* `testInterProGeneOntology.py` to determine which InterPro domains imply a given GO term (plus descendants).
- **`go-basic.obo`** — GO ontology in OBO format.  
  *Used by:* `testInterProGeneOntology.py` to compute **`is_a` descendants** of the target GO term.

> The GO module will emit friendly `[info]/[warn]` messages if any of the files are missing. With an empty/missing map or DAG, matches will tend to be `False` (as expected).


---

## 4) Core modules — what they do and how to run them

### 4.1 `queryWbpBiomart.py` — central client
- **Purpose:** single function that sends an XML to WBP BioMart and returns a DataFrame.
- **Functions:**
  - `fetch_wbp_biomart_using_xml_polars(xml_query)` → `pl.DataFrame`
  - `fetch_biomart_with_validation(xml_query, expected_columns)` → `pl.DataFrame` (validates headers immediately)
- **Behaviour:**
  - Retries with exponential backoff.
  - Enforces non-empty responses.
  - Parses as **Polars** by default.

**No direct user invocation** — it’s imported by the other modules.

---

### 4.2 `populateGeneList.py` — principal gene list (per species)
- **Purpose:** Get the **protein-coding** gene list for a single species and select a **principal isoform** per gene (longest CDS, with deterministic ties).
- **XML:** headered (`header="1"`) and uses attributes for: gene ID/name, transcript ID/biotype, description, CDS start/end.  
- **Key functions:**
  - `retrieveGeneListFromWbpBiomart(species_code) -> pd.DataFrame`  
    - Validates **single species**.
    - Fetches **genes** and **CDS** tables via `queryWbpBiomart.py` with **explicit expected columns** checks.
    - Computes CDS length and **selects principal isoform** deterministically.
    - **Returns Pandas** with BioMart labels preserved.
- **Output columns (subset):**
  - `"Gene stable ID"`, `"Transcript stable ID"`, `"Genome name"`, `"Gene description"`, `"Transcript biotype"`
  - `"principal_transcript_id"`, `"principal_cds_length"`, `"principal_isoform_type"`
  - plus a `"species"` column for bookkeeping.
- **Integrity:** asserts no duplicate `"Gene stable ID"` in the per-species result.

**Example (from repository root):**
```bash
python -c "from populateGeneList import retrieveGeneListFromWbpBiomart as run; df = run('trtricprjeb535'); print(df.head().to_string()); print(len(df))"
```

---

### 4.3 `testInterProGeneOntology.py` — InterPro → GO test
- **Purpose:** Flag genes whose InterPro domains imply a **target GO term** or any of its **`is_a` descendants** (via `interpro2go` + `go-basic.obo`). Returns a boolean `<test_name>` plus evidence text.
- **Contract details:**
  - **Exact XML** from the original, with only `header="1"` and `{species_code}` interpolated.
  - Validates **single species**.
  - **Header discovery prints** the discovered column labels.
  - **No renaming** of BioMart labels; expected columns are:
    - `"Gene stable ID"`, `"InterPro ID"`, `"InterPro description"`
  - Computes descendants with **`is_a` only** (matches original behaviour).
  - **Output:** `"Gene stable ID"`, `<test_name>` (bool), `<test_name>_evidence` (string).
  - **Integrity:** asserts **no duplicate** `"Gene stable ID"` in its output.
- **Required local files:** `interpro2go`, `go-basic.obo` (same folder).
- **Example:**
```bash
python -c "from testInterProGeneOntology import testInterProGeneOntology as run; df = run('trtricprjeb535','GO:0003824','is_enzyme'); print(df.head(10).to_string())"
```

---

### 4.4 `testInterProStringSearch.py` — InterPro regex test
- **Purpose:** Flag genes whose InterPro annotations **text-match** a given regex (e.g., GPCR). Returns `<test_name>` boolean + evidence text.
- **Contract details:**
  - **Exact XML** (only `header="1"` + `{species_code}`).
  - Validates **single species**.
  - **Header discovery prints** the labels.
  - Required columns: `"Gene stable ID"`, `"InterPro ID"`, `"InterPro description"`.
  - Builds `"interpro_annotation" = "InterPro ID - InterPro description"` and matches a **Polars regex**.
  - **Output:** `"Gene stable ID"`, `<test_name>` (bool), `<test_name>_evidence`.
  - **Integrity:** asserts **no duplicate** `"Gene stable ID"` in its output.
- **Example (GPCR, original pattern):**
```bash
python -c "from testInterProStringSearch import testInterProStringSearch as run; df = run('trtricprjeb535', r'G-protein coupled receptor|GPCR.[^k]', 'is_gpcr'); print(df.head(10).to_string())"
```
> You can widen the regex **argument** later if required (e.g., case-insensitive variants), without touching module code.


---

## 5) Orchestrator — `run_v01.py`

`run_v01.py` is the **only entry point** you need to run. It creates a timestamped run folder, writes logs, records metadata, calls all modules per species, joins outputs, and writes `export.csv`/`export.json`.

### 5.1 Species list
Open `run_v01.py` and edit:
```python
SPECIES_LIST = ["trtricprjeb535", "wubancprjna275548"]  # example can run one or multiple
```

### 5.2 What it does
1. Creates `pipeline_runs/<timestamp>_<LABEL>/` and starts a tee logger (`log.txt`) and `run.json` with metadata.
2. For each species:
   - Calls **principal gene list** (`populateGeneList.retrieveGeneListFromWbpBiomart`).
   - Calls **human orthologues** module (keeps BioMart labels; merged left on `"Gene stable ID"`).
   - Calls **InterPro GO tests**: e.g., `is_enzyme` (`GO:0003824`), `is_kinase` (`GO:0004672`).
   - Calls **InterPro string test**: e.g., `is_gpcr` with the **original** pattern `G-protein coupled receptor|GPCR.[^k]`.
   - Left‑merges each module output on `"Gene stable ID"`; booleans filled to `False` if missing.
3. Concatenates per‑species tables into one `combined` table.
4. Updates `run.json` counts (`total_rows`, per‑species counts).
5. Writes:
   - `export.csv` (full combined table)
   - `export.json` (records array)

### 5.3 Outputs
Inside `pipeline_runs/<timestamp>_<LABEL>/` you’ll see:
- `log.txt` — full console log (includes **header discovery** prints from modules).
- `run.json` — run metadata and per‑species row counts.
- `export.csv`, `export.json` — final combined export.

### 5.4 Sanity asserts (recommended)
Add these tiny checks (already suggested) right before writing `export.csv/json`:
- **No duplicate** `"Gene stable ID"` in the final table.
- Ensure any present flags (`is_enzyme`, `is_kinase`, `is_gpcr`) are **booleans** and each has a corresponding `*_evidence` column.
These are small `assert`/`if` blocks and print final TRUE counts to the log.


---

## 6) Running the whole pipeline

> **Prereqs:** Place `interpro2go` and `go-basic.obo` in the same folder as the modules.

### Windows PowerShell
```powershell
cd path\to\pipeline
python run_v01.py
```

### macOS / Linux
```bash
cd /path/to/pipeline
python3 run_v01.py
```

**What you should see:**
- Early lines with run folder and species list.
- Per-module header discovery like:
  ```
  [info] Biomart header discovered (3 columns):
    [01] Gene stable ID
    [02] InterPro ID
    [03] InterPro description
  ```
- Per-species row counts and final `Combined export shape`.
- Final `export.csv`/`export.json` and `run.json` in the run folder.


---

## 7) Data schema — common columns

Exact columns depend on which modules you wire in, but the **canonical keys/labels** are:

**Principal gene list (`populateGeneList.py`):**
- `"Gene stable ID"` (key), `"Transcript stable ID"`, `"Genome name"`, `"Gene description"`, `"Transcript biotype"`
- `"principal_transcript_id"`, `"principal_cds_length"`, `"principal_isoform_type"`

**Human orthologues module (wired in `run_v01.py`):**
- `"Human gene stable ID"`, `"Human gene name"`, `"Homology type"`
- `"% identity"`, `"Human % identity"`
- `"lacks_WBP_human_orthologue"`, `"lacks_WBP_human_orthologue_evidence"`
- `"best_WBP_human_orthologue_lt_40pct_identity"`, `"best_WBP_human_orthologue_lt_40pct_identity_evidence"`

**InterPro GO module (`testInterProGeneOntology.py`):**
- `"Gene stable ID"`, `<test_name>`, `<test_name>_evidence`  
  e.g., `is_enzyme` / `is_enzyme_evidence`, `is_kinase` / `is_kinase_evidence`.

**InterPro string module (`testInterProStringSearch.py`):**
- `"Gene stable ID"`, `<test_name>`, `<test_name>_evidence`  
  e.g., `is_gpcr` / `is_gpcr_evidence`.

> All columns use the **exact BioMart labels** as discovered in the headered response.


---

## 8) Reproducing legacy exports (parity tips)

- Use the **original GPCR regex** first: `G-protein coupled receptor|GPCR.[^k]`.  
  If you later need broader coverage (e.g., case‑insensitive, “7TM”, “serpentine”), change **only the `search_string` argument** in `run_v01.py`.
- Confirm your **species list** matches the legacy export (counts add up across species).  
- Keep `interpro2go` and `go-basic.obo` versioned in your repo to avoid silent drift.


---

## 9) Troubleshooting

- **“Query ERROR: …” appears as the only column**  
  Your XML was rejected by BioMart. Check that the template matches the original, **only** changing `header="1"` and `{species_code}`.

- **Missing expected columns**  
  The modules will fail early and print the actual header list — adapt the XML attributes if needed (after explicit approval).

- **All `False` for GO module**  
  Ensure `interpro2go` and `go-basic.obo` exist in the same folder. The module logs how many links/nodes were loaded and how many descendants were considered.

- **Duplicates detected**  
  The final sanity asserts raise an error mentioning sample `"Gene stable ID"` values. Investigate earlier merges.


---

## 10) Extending the pipeline — writing new modules

When adding new analysis modules, follow this minimal scaffold:

1. **Exact XML** (headered, single species) — only `header="1"` and `{species_code}` substitutions are allowed.
2. **Fetch with central client** → immediately **print header discovery**.
3. **Validate minimal required columns** for your logic.
4. **Polars internally**; **return Pandas** with columns:
   - `"Gene stable ID"`, `<test_name>`, `<test_name>_evidence`
5. **Determinism & integrity** — assert no duplicate `"Gene stable ID"` before returning.
6. **Join in `run_v01.py`** with a left-merge on `"Gene stable ID"`, fill booleans to `False`.

This keeps the codebases consistent and easy to reason about.


---

## 11) Quick examples

**Single-species GPCR table (ad‑hoc):**
```bash
python -c "from testInterProStringSearch import testInterProStringSearch as run; df = run('trtricprjeb535', r'G-protein coupled receptor|GPCR.[^k]', 'is_gpcr'); print(df['is_gpcr'].sum())"
```

**Single-species enzyme table:**
```bash
python -c "from testInterProGeneOntology import testInterProGeneOntology as run; df = run('trtricprjeb535','GO:0003824','is_enzyme'); print(df.head(10).to_string())"
```

**Full pipeline run:**
```bash
python run_v01.py
```


---

## 12) FAQ

- **Where do I put `interpro2go` and `go-basic.obo`?**  
  In the **same folder** as the modules (`pipeline/`). The GO module looks **only** in its own directory.

- **Do I need goatools?**  
  No. The GO module performs a **minimal OBO parse** (`is_a` edges only).

- **Can I run multiple species at once?**  
  Yes—edit `SPECIES_LIST` in `run_v01.py`. Each species is processed independently and concatenated.

- **Why are there fewer rows in module outputs than the final export?**  
  Module outputs include only genes with InterPro rows. The orchestrator left‑joins module outputs onto the full principal list, so the final export has all genes; flags default to `False` when no InterPro data exists.


---

## 13) License / provenance

- Module internals follow your explicit contract and keep all **BioMart labels verbatim**.
- `run_v01.py` owns run folder creation, logging (`log.txt`), and `run.json` with per‑species counts.
- Keep versioned copies of `interpro2go` and `go-basic.obo` to ensure reproducibility across runs.
