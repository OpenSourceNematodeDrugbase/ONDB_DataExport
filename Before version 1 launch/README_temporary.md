
# Open-Source Nematode DrugBase (ONDB) — Full Pipeline Documentation 

## Overview
The Open-Source Nematode DrugBase (ONDB) is a fully reproducible, contract‑driven data pipeline designed to integrate genomic, functional, orthology, and phenotypic information across parasitic nematodes. The system unifies **BioMart-derived structural and orthology data** with **WormMine-derived RNAi phenotypic evidence** from *Caenorhabditis elegans* to support target prioritisation for neglected tropical disease drug discovery.

This README documents the full architecture as of **23 November 2025**, including new modules, data flow, contracts, provenance, and the integration strategy for C. elegans RNAi phenotypes.

---

# 1. Architectural Structure

The pipeline is divided into three major subsystems:

## **1. BioMart modules (structural + orthology + protein domains)**  
These fetch all genomic and functional annotations required for each parasite species:
- Principal gene list  
- Coding sequence (CDS) boundaries  
- Human orthologues  
- *C. elegans* orthologues  
- InterPro domains  
- InterPro→GO mapping  
- GO DAG propagation  

## **2. WormMine modules (functional phenotype extraction)**  
This subsystem retrieves RNAi phenotypes from curated *C. elegans* RNAi screens via the WormMine REST interface:
- Uses template-based XML queries identical to the WormBase WormMine QueryBuilder  
- Produces long-format RNAi tables  
- Aggregates phenotypes into per-gene flags (e.g., lethal, sterile, locomotion)

## **3. Integration subsystem (run_v01.py)**  
This orchestrates:
- Species‑wise execution of all BioMart and WormMine modules  
- Validation and deterministic merging  
- Per-gene consolidation  
- Final export table  
- Provenance metadata (run folders, URLs, version numbers, logs)

The following diagram summarises the data flow:

```
                    +------------------+
                    |  Species list    |
                    +------------------+
                              |
                --------------------------------
                |                              |
         +-------------+                +----------------+
         |   BioMart   |                |   WormMine     |
         |  (per-spec) |                | (C. elegans)   |
         +-------------+                +----------------+
         | Gene list   |                | RNAi XML       |
         | CDS         |                | REST query     |
         | Human orth. |                | Long-format    |
         | CE orth.    |                | Aggregation    |
         | InterPro    |                +----------------+
         | GO mapping  |                        |
         +-------------+                        |
                |                                |
                ------------ merged ---------------
                              |
                     +---------------------+
                     |     Final export    |
                     |    (per species)    |
                     +---------------------+
```

---

# 2. Core Development Contract (ONDB Pipeline Contract)

All ONDB modules must satisfy the following rules:

### **2.1 BioMart XML consistency**
- XML templates must be used exactly as provided by BioMart.
- The only allowed modifications:
  1. `header="1"` (mandatory)
  2. Insert `{species_code}` into species filter attributes
- Do NOT modify any attribute order, whitespace, or structure.

### **2.2 Column discovery and validation**
- Expected column names must come **from the XML template itself**, not hard-coded.
- After fetching, the module must:
  - Print the header line  
  - Print a numbered list of column names  
  - Validate required columns with alias‑based drift correction  

### **2.3 Deterministic row selection**
- One row per gene after merging.  
- Ties resolved deterministically (highest % identity, then lexical fallback).

### **2.4 Single-species input**
- Species codes must be supplied exactly one at a time.
- Comma‑separated species codes are rejected.

### **2.5 Polars internal → Pandas boundary**
- All internal processing uses **Polars** for performance.
- Each module returns **Pandas** to the boundary for merging.

### **2.6 Provenance**
- Downloaded auxiliary files (InterPro→GO, GO DAG) are stored in the run folder.
- Versions and URLs logged into `run.json`.

---

# 3. BioMart Modules

Below is the updated full documentation matching today’s pipeline state.

---

## **3.1 Principal Gene List**

### Purpose  
Retrieve major structural gene information per species.

### Source  
BioMart: `wbps_gene` dataset.

### XML Attributes  
- `Gene stable ID`
- `Transcript stable ID`
- `Transcript biotype`
- `Gene description`
- `Genome name`

### Validation  
- Required columns discovered at runtime.  
- Only protein-coding transcripts retained.  
- Principal isoforms selected via deterministic rule:
  - Longest CDS
  - If equal, pick lowest transcript ID lexically

### Output  
One row per parasite gene.

---

## **3.2 Coding Sequences (CDS)**

### Purpose  
Retrieve CDS boundaries for defining principal isoforms.

### Attributes  
- `CDS start (within cDNA)`  
- `CDS end (within cDNA)`  

### Rules  
- Multi‑exonic CDS reconstructed from Polars groupby.  
- Deterministic isoform selection.

---

## **3.3 Human Orthologues**

### XML attributes  
- `Human gene stable ID`  
- `Human gene name`  
- `Homology type`  
- `% identity`  
- `Human % identity`  

### Rules  
- Aliased columns allowed (`Human % identity` sometimes becomes `% identity r1`).  
- One row per parasite gene post‑resolution.

---

## **3.4 C. elegans Orthologues**

### XML attributes  
- `Caenorhabditis elegans (PRJNA13758) [WS290] gene stable ID`  
- `gene name`  
- `Homology type`  
- `% identity`  
- CE‑specific % identity  

### Output  
Added to export as:
- `ce_gene_id`  
- `ce_gene_name`  
- `best_ce_identity_homology_type`  
- Flags:
  - `lacks_celegans_orthologue`
  - `best_celegans_orthologue_lt_40pct_identity`

---

## **3.5 InterPro Domains**

### XML attributes  
- `InterPro ID`  
- `InterPro description`  

Multiple rows per gene allowed.

---

## **3.6 GO Annotations**

### Inputs  
- `interpro2go` mapping (IPR → GO)  
- GO DAG (is_a hierarchy)

### Outputs  
Columns with boolean flags:
- `is_enzyme`  
- `is_kinase`  
- `is_gpcr`  

---

# 4. WormMine Modules (NEW, 23 Nov 2025)

## **4.1 Motivation**

*C. elegans* provides the most comprehensive genome‑wide RNAi phenotyping available. WormMine exposes these data via a declarative XML query format. Integrating these phenotypes allows ONDB to prioritise parasite genes supported by functional evidence.

### Functional questions supported:
- Is the orthologue essential?  
- Does knockdown cause sterility, embryonic lethality, locomotion defects?  
- Are phenotypes consistent across multiple RNAi screens?  

---

## **4.2 XML Template for WormMine (canonical)**

ONDB uses *the exact XML query downloaded from WormMine*:

```xml
<query model="genomic" view="
    Gene.primaryIdentifier
    Gene.secondaryIdentifier
    Gene.RNAiResult.primaryIdentifier
    Gene.RNAiResult.phenotype.identifier
    Gene.RNAiResult.phenotype.name
    Gene.RNAiResult.phenotypeRemark
    Gene.RNAiResult.remark"
    sortOrder="Gene.primaryIdentifier ASC">
  <constraint path="Gene.primaryIdentifier" op="=" value="WBGene00006757" code="A" />
</query>
```

Only the value of `Gene.primaryIdentifier` is changed programmatically.

---

## **4.3 WormMine REST API (Python 3.12–compatible)**

Because the official `intermine` Python package is Python 2–oriented and incompatible with Python 3.12 (`urlparse`, `MutableMapping` deprecations), ONDB uses a **custom REST client**:

- Sends XML payload to:
  ```
  http://intermine.wormbase.org/tools/wormmine/service/query/results?format=json
  ```
- Robust against:
  - HTTP 500 errors
  - Timeouts
  - Restart of server connection
  - Empty results

---

## **4.4 Long‑Format Results**

Each RNAi result is returned as one row:

Columns:
- `wormbase_gene_id`
- `sequence_name`
- `rnai_result_id`
- `phenotype_id`
- `phenotype_name`
- `phenotype_remark`
- `rnai_result_remark`

Multiple rows per CE gene are expected.

---

## **4.5 Batch Parallel Strategy C (Default for ONDB)**

### Motivation  
WormMine handles many small requests better than one huge request.

### Strategy  
1. Deduplicate CE gene IDs  
2. Break into batches of size 20  
3. Sequential batches, threads inside each batch  
4. Truncated exponential backoff for retries  

---

## **4.6 Aggregation per C. elegans Gene**

The long table is collapsed to one row per CE gene.

### Phenotypic flags generated:
- `ce_rnai_any_lethal`  
- `ce_rnai_any_emb` (embryonic lethal)  
- `ce_rnai_any_sterile`  
- `ce_rnai_any_locomotion`  

Rules:
- A phenotype is considered present if **any** RNAi experiment mapped to that category.

Unused phenotype names are kept in a parallel list column for future search functionality.

---

## **4.7 Missing CE orthologues**

If a parasite gene lacks valid CE orthologues:
- No WormMine call is made  
- All CE RNAi flags default to `False`  

---

# 5. Integration Subsystem (run_v01.py)

## Full process for each species

### **5.1 Stepwise workflow**
1. Load species list  
2. For each species:
   - Gene list  
   - CDS  
   - Human orthologues  
   - CE orthologues  
   - InterPro / GO  
   - Merge all  
   - CE RNAi selection:
     - Extract CE gene IDs  
     - Fetch RNAi data  
     - Aggregate  
   - Attach phenotype flags  
3. Concatenate all species  
4. Export final table + run folder metadata  

---

## 5.2 Error Handling
- Column drift warnings  
- BioMart unavailable  
- XML mismatch  
- WormMine 500 errors  
- Timeout → retry  
- Missing orthologues  

---

## 5.3 Provenance Logging
Each run folder contains:
- `run.json`  
- raw BioMart XML queries  
- exact URLs  
- InterPro→GO copies  
- GO DAG copy  
- outputs per species  
- combined export  

---

# 6. Final Export Format

### Columns (example, 42 columns)
- Structural columns  
- Human orthology columns  
- CE orthology columns  
- CE RNAi phenotype flags  
- Functional domain flags  

Each row = one parasite gene.

---

# 7. Known Limits and Future Work
- WormMine rate limiting still variable  
- Phenotype grouping: future expansion to full ontology  
 

---

# 8. Citation and Credits
The ONDB pipeline is designed for transparent and reproducible target discovery, integrating open data from:
- WormBase ParaSite BioMart  
- WormBase WormMine  
- Gene Ontology  
- InterPro  
- ZINC/ChEMBL future integration  

