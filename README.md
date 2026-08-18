# PeptiLine: An Interactive Platform for Customizable Functional Peptidomic Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Django 4.2](https://img.shields.io/badge/django-4.2-092E20.svg)](https://www.djangoproject.com/)
[![DOI](https://img.shields.io/badge/DOI-ZENODO__DOI__PLACEHOLDER-blue.svg)](https://doi.org/ZENODO_DOI_PLACEHOLDER)

PeptiLine turns peptidomic mass spectrometry output into annotated, statistically compared, and sequence-mapped visualizations, without requiring installation, a programming environment, or a local copy of any bioactivity database.

- **Live application:** https://mbpdb.nws.oregonstate.edu/peptiline/ (no account or install required)
- **Publication:** *Journal of Proteome Research* (ACS), MS pr-2025-01102w, under review
- **License:** MIT
- **Archived version:** Zenodo DOI `ZENODO_DOI_PLACEHOLDER` (inserted on code deposit; see [Data and reproducibility](#data-and-reproducibility))

## Table of Contents
- [Overview](#overview)
- [Two ways to use PeptiLine](#two-ways-to-use-peptiline)
- [Modules](#modules)
  - [1. Data Transformation](#1-data-transformation)
  - [2. Data Analysis](#2-data-analysis)
  - [3. Heatmap Visualization](#3-heatmap-visualization)
- [Input data requirements](#input-data-requirements)
- [Installation (local development)](#installation-local-development)
- [MBPDB integration](#mbpdb-integration)
- [Data and reproducibility](#data-and-reproducibility)
- [System and dependency versions](#system-and-dependency-versions)
- [Project structure](#project-structure)
- [Documentation](#documentation)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)
- [Authors](#authors)
- [Contributing](#contributing)
- [Legacy notebook implementation](#legacy-notebook-implementation)

## Overview

PeptiLine is a Django web application with three linked dashboards: **Data Transformation**, **Data Analysis**, and **Heatmap Visualization**. A dataset prepared in Data Transformation carries directly into either visualization dashboard, so a full analysis runs as one continuous session rather than three disconnected tools. Background processing (Celery + Redis) handles longer transformations and BLAST searches without blocking the browser.

The hosted deployment queries the Milk Bioactive Peptide Database (MBPDB) directly, so functional annotation works with no local database setup. The application code in this repository is independent of MBPDB's database and search backend (see [MBPDB integration](#mbpdb-integration)) and can be run standalone against a user-supplied functional annotation table.

## Two ways to use PeptiLine

**Hosted (recommended):** https://mbpdb.nws.oregonstate.edu/peptiline/. Full functionality including live MBPDB search, no installation.

**Local:** clone this repository and run the Django app yourself (see [Installation](#installation-local-development)). Everything works locally except live MBPDB search, which depends on a private database this repository does not include; upload a pre-downloaded MBPDB TSV or your own functional annotation table instead (see [MBPDB integration](#mbpdb-integration)).

## Modules

### 1. Data Transformation

Ingests raw or pre-annotated peptidomic data and prepares it for analysis.

- Multi-format upload (CSV, TSV, TXT, XLSX) with automatic column mapping for Proteome Discoverer, MaxQuant, Skyline, PEAKS, Spectronaut, and native PepEx exports; inline modifications (`PEP(+57.02)TIDE`, `PEPTIDE[DN]`, enzymatic flanks) are stripped so sequences stay valid for exact/homology matching.
- Automatic long-format detection for engines that export one row per PSM or precursor (PEAKS, Spectronaut): rows are pivoted to a peptide-by-sample matrix and summed/collapsed within a sample before further processing.
- Multi-dataset import: upload several files of the same format at once; PeptiLine merges them on shared peptide identity, keeps each file's abundance columns, and pads missing peptides with `NA`.
- Functional annotation either by MBPDB search (exact match or a user-set homology threshold) or by uploading a custom function table (see [Table 2](#table-2-functional-annotation-columns)).
- Multi-protein peptide resolution: split, retain, or remove ambiguous multi-protein assignments; a **Merge/Rename Protein Sources** control folds selected accessions, including single proteins, into one canonical name (for example, consolidating separately searched β-casein A1/A2 into one protein).
- Study-variable / group assignment via the UI or a JSON upload, including hierarchical grouping (e.g., a "Bitter" group composed of two finer-grained groups).
- Reusable protein-mapping and column-rename keys: mapping and renaming decisions can be exported as JSON and re-applied to new datasets without repeating the manual steps.
- Twelve export types, including the merged dataset, MBPDB search results, summed functional data, sample-to-sample and biological-replicate correlation tables (Pearson/Spearman, optional log10), per-group peptide lists, and the two reusable keys above.

### 2. Data Analysis

Interactive descriptive and comparative statistics on the transformed dataset.

- Grouped bar, stacked bar, pie, and correlation scatter-matrix plots, viewable by sample, protein, or function, in absolute or relative terms, for peptide count or abundance.
- Replicate-based SEM error bars on grouped bar plots.
- Significance testing on grouped bar plots (absolute values): one-way ANOVA with Tukey HSD by default, or Welch's ANOVA with Games-Howell for heteroscedastic data. Two-group comparisons render as a bracket annotation; three or more groups render as compact letter displays. Requires at least three replicates per group.
- Filtering by protein, bioactive function, or both, with an option to separate functional from non-functional peptides.
- Export as interactive HTML, or as static PNG/SVG rendered server-side at 600 dpi, matching what is shown on screen.

### 3. Heatmap Visualization

Maps peptide density and abundance onto a protein's amino-acid sequence.

- Landscape, portrait, and compact orientations from a single rendering engine, so all three stay behaviorally consistent.
- A differential comparison track overlays a signed effect size, Cohen's *d* or log2 fold change, between two sample groups on one protein or between two protein variants within one group, with small/medium/large effect-size reference thresholds. Requires replicate-level input.
- Optional amino-acid lettering printed directly on each sample's tile, for residue-level comparison across proteins or variants.
- Strip-start-sequence option that trims a protein's N-terminal signal peptide, either from UniProt annotation or a manual residue-count override, so the map reflects the mature protein.
- UniProt lookup for reference sequences, or upload a custom FASTA for protein variants not distinguished in UniProt (for example, β-casein A1 vs. A2).
- Export as interactive HTML, or as static PNG/SVG rendered server-side at 600 dpi.

## Input data requirements

### Table 1. Required columns for peptidomic data

| Data type | Name of created column | Acceptable source columns |
|-----------|------------------------|-------------------|
| Peptide sequence | Sequence, Unique Peptide ID | peptide, sequence, Annotated Sequence, Peptide Sequence |
| Precursor protein ID | Protein | Leading razor protein, UniProt ID, protein, Proteins, Protein, prot_acc, Accession, Master Protein Accessions, Protein ID |
| Peptide start | start | start position, start, Start, pep_res_before, Positions in Master Proteins |
| Peptide end | end | end position, end, End, pep_res_after, Positions in Master Proteins |
| Sample intensity | user selected | user selected |
| Modifications | Unique Peptide ID | Modified Sequence, Modifications, modified_peptide |

Native exports from Proteome Discoverer, MaxQuant, Skyline, PEAKS, and Spectronaut are auto-detected without hand-editing. FragPipe columns are named in the mapping table but not yet verified against a real export; treat FragPipe support as unconfirmed until an export has been tested.

### Table 2. Functional annotation columns

| Required column | Explanation |
|------------------|-------------|
| search_peptide | Matches your peptidomic dataset's sequence exactly; the key linking functional and peptidomic data. |
| peptide | The database sequence matched, which may differ slightly under homology matching. |
| function | The bioactivity label used for grouping and annotation. |

### Example files

Eight example files ship with the repository under `examples/`: a peptidomic dataset, a merged/transformed dataset, an MBPDB-format functional annotation table, a protein FASTA, and one example each of the study-variables, column-rename, protein-mapping, and technical-replicate JSON keys. These are illustrative, not case-study data; see [Data and reproducibility](#data-and-reproducibility) for the manuscript's own dataset.

## Installation (local development)

```bash
# Clone the repository
git clone https://github.com/kuhfeldrf/peptiline.git
cd peptiline

# Create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt

# Apply migrations and run
python manage.py migrate
python manage.py runserver
```

Background tasks (large transformations, exports) require a running Celery worker and Redis broker; see `docs/INSTALL.md` for the full local and Docker setup, including nginx configuration for a production-style deployment.

MBPDB search will not run locally (see [MBPDB integration](#mbpdb-integration)); every other feature in all three modules is fully functional against local data.

## MBPDB integration

**MBPDB database search is not included in this repository.** MBPDB is a separate, privately maintained database; the hosted PeptiLine deployment queries it directly, but the search backend and underlying records are not part of this codebase.

**Options for local use:**
1. Use the hosted app for MBPDB search, then download the results as TSV and upload them here for further analysis.
2. Upload your own functional annotation table (see [Table 2](#table-2-functional-annotation-columns)).
3. Use the bundled `examples/example_mbpdb_functional_data.csv` to try the workflow without a real search.

Everything else, including data transformation, UniProt sequence retrieval, statistics, and both visualization modules, runs identically locally and on the hosted deployment.

## Data and reproducibility

The manuscript's case study (bitterness in aged Cheddar cheese) used PeptiLine to generate Figures 5-7 and Supplementary Figures S1-S6. Raw mass spectrometry data are deposited at ProteomeXchange, accession **PXD079655**. The exact application version used to generate the manuscript's figures is archived on Zenodo: `ZENODO_DOI_PLACEHOLDER`.

Supporting tables and figures (S1-S13, S1-S6) are provided in `supplementals/`, matching what the hosted app serves at `/peptiline/supplementals/`. A step-by-step reproduction protocol, module-by-module settings for every case-study figure, is in `docs/REPRODUCIBILITY.md`.

## System and dependency versions

- **Python:** 3.10+
- **Framework:** Django 4.2
- **Background processing:** Celery, Redis
- **Plotting/export:** Plotly, Kaleido (server-side static export, bundled headless Chromium, no system Chrome dependency), Pillow
- **Statistics:** SciPy, statsmodels
- **Data handling:** pandas, NumPy
- **Sequence search:** NCBI BLAST+ (`blastp`), required only for the MBPDB-search code path
- **OS:** Linux, macOS, or Windows via Docker; the hosted deployment runs on Ubuntu Linux

Full pinned versions in `requirements.txt`.

## Project structure

```
peptiline/
├── data_transformation/   # Module 1: upload, mapping, annotation, export
├── data_analysis/         # Module 2: descriptive statistics and plots
├── heatmap_viz/            # Module 3: sequence-mapped heatmaps
├── utils/                 # Shared helpers (UniProt client, static export)
├── examples/               # Sample input files
├── supplementals/          # Manuscript Supporting Information (tables + figures)
├── docs/                   # USER_GUIDE, REPRODUCIBILITY, TROUBLESHOOTING, INSTALL
├── requirements.txt
├── LICENSE
└── README.md               # This file
```

## Documentation

- **Feature and workflow guide:** `docs/USER_GUIDE.md`
- **Case-study reproduction protocol:** `docs/REPRODUCIBILITY.md`
- **Local install and deployment:** `docs/INSTALL.md`
- **Troubleshooting:** `docs/TROUBLESHOOTING.md` (or the summary below)

## Troubleshooting

**Timeouts or script errors in the browser.** Refresh the page; this resolves most transient timeout or script-execution issues.

**Memory errors on large datasets.** Close other applications, process data in smaller batches, or use a machine with more RAM for datasets above roughly 10,000 peptides.

**UniProt connection timeouts.** Check your internet connection, or use local files instead: upload a FASTA directly rather than relying on the live UniProt lookup.

**Python version errors.** Confirm `python --version` reports 3.10 or higher, and that the virtual environment is active before reinstalling dependencies.

## Citation

Manuscript under review at the *Journal of Proteome Research* (ACS); citation will be updated on acceptance.

```
Kuhfeld R, Nielsen SD-H, Dallas DC. PeptiLine: An Interactive Platform for
Customizable Functional Peptidomic Analysis. Submitted to the Journal of Proteome Research. 2026.
```

**BibTeX:**
```bibtex
@article{kuhfeld2026peptiline,
  title={PeptiLine: An Interactive Platform for Customizable Functional Peptidomic Analysis},
  author={Kuhfeld, Russell and Nielsen, S{\o}ren D-H and Dallas, David C},
  journal={Journal of Proteome Research},
  year={2026},
  note={Submitted, under review}
}
```

## License

MIT License. See [LICENSE](LICENSE).

## Authors

- Russell Kuhfeld (Oregon State University)
- Søren D-H. Nielsen (Arla Foods Ingredients Group P/S)
- David C. Dallas (Oregon State University)

## Contributing

Contributions are welcome:
1. Fork the repository
2. Create a feature branch
3. Open a pull request describing the change

For larger changes, open an issue first to discuss scope.

## Legacy notebook implementation

The original Jupyter/ipywidgets notebook implementation (v1, accompanying the initial manuscript submission) is archived at `github.com/kuhfeldrf/peptiline`, tagged `v1.0`. It is kept for historical reference and is not maintained; use this repository for current development and for the deployed application.
