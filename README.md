# PeptiLine: An Interactive Platform for Customizable Functional Peptidomic Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)

> **⚠️ Legacy Repository Notice**
>
> This repository contains the **original Jupyter notebook implementation** of PeptiLine (v1.0).
> The pipeline has since been converted into a fully deployed Django web application with an
> improved UI, background task processing, and direct MBPDB database integration.
>
> **Production web app source code:** https://github.com/mbpdb/mbpdb
> **Live application:** https://mbpdb.nws.oregonstate.edu/peptiline/
>
> These notebooks remain available for local use, reproducibility of the original manuscript
> figures, and as a reference implementation. See the notes below regarding deprecated widgets.

## Overview
PeptiLine transforms peptidomic mass spectrometry data into interactive visualizations and statistical summaries.

> **Note:** The interactive Jupyter notebooks (`heatmap_visualization.ipynb`, `data_analysis.ipynb`) that used ipywidgets are no longer supported. The web application at https://mbpdb.nws.oregonstate.edu/peptiline/ provides equivalent functionality via the Data Transformation, Data Analysis, and Heatmap dashboards.

- **Publication:** Submitted to PLOS Computational Biology (2025, under review)
- **Web Application:** https://mbpdb.nws.oregonstate.edu/peptiline/
- **License:** MIT (OSI compliant)
- **Note:** MBPDB database search is only available via web application. See [Important Note](#important-note-mbpdb-search-functionality).

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Testing the Installation](#testing-the-installation)
- [Getting Started](#getting-started)
  - [Input Data Requirements](#input-data-requirements)
- [Reproducibility](#reproducibility)
- [External Dependencies](#external-dependencies)
- [Documentation](#documentation)
- [Dependencies](#dependencies)
- [Project Structure](#project-structure)
- [Support and Contact](#support-and-contact)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)
- [Authors](#authors)
- [Version History](#version-history)
- [Contributing](#contributing)
- [Important Note: MBPDB Search Functionality](#important-note-mbpdb-search-functionality)
- [Software Archive](#software-archive)

## Features
1. **Data Transformation** - Organizes and annotates peptidomic data with multi-protein peptide assignment and MBPDB bioactivity annotation
2. **Descriptive Analysis** - Interactive visualizations including correlation plots, functional distributions, and protein origin contributions
3. **Heatmap Visualization** - Peptide density and positioning maps along protein sequences

## System Requirements
- **OS:** Windows 10+, macOS 11+, Linux (Ubuntu 20.04+)
- **Python:** 3.12+
- **Jupyter:** 7.4.5+
- **RAM:** 4 GB minimum (8 GB for >5000 peptides)
- **Storage:** 500 MB + user data
- **Internet:** Required for installation, optional for UniProt API

## Installation

### Web-Based Access
https://mbpdb.nws.oregonstate.edu/peptiline/ (includes MBPDB database search)

### Local Installation
See Quick Start in [Getting Started](#getting-started). Note: MBPDB search unavailable locally (see [Important Note](#important-note-mbpdb-search-functionality)).

## Testing the Installation

Verify installation using provided example data (5-10 minutes total):

1. **Launch Jupyter Lab:**
   ```bash
   cd peptiline
   source .venv/bin/activate  # Windows: .venv\Scripts\activate
   jupyter lab
   ```

2. **Test each module:**
   - Open `data_transformation.ipynb`, upload `examples/example_peptide_data.csv`, run cells
   - Open `data_analysis.ipynb`, upload transformed dataset, run cells  
   - Open `heatmap_visualization.ipynb`, select any protein, generate heatmap

**Success:** All cells execute without errors, widgets display, visualizations render.

See [Troubleshooting](#troubleshooting) if issues occur.

## Getting Started

### Using the Notebooks Locally
1. After launching Jupyter Lab, navigate to the desired module:
   - `data_transformation.ipynb` - Start here for data processing
   - `data_analysis.ipynb` - For descriptive analysis and visualizations
   - `heatmap_visualization.ipynb` - For sequence-based heatmap generation

2. Run the first cells to initialize the interactive interface

3. Follow the on-screen instructions within each notebook

### Quick Start (Recommended)
```bash
# Clone the repository
git clone https://github.com/Kuhfeldrf/peptiline.git
cd peptiline

# Create and activate a virtual environment (recommended)
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install --upgrade pip
pip install -r notebook_requirements.txt

# Launch Jupyter Lab
jupyter lab
```

### Input Data Requirements
PeptiLine accepts common peptidomic file formats:
- **Peptidomic data:** CSV, TSV, TXT, or XLSX files from Proteome Discoverer, MaxQuant, Skyline, PEPEX

   - Table 1. Required Columns for Peptidomic Data

      | Data type | Name of created columns | Acceptable columns |
      |-----------|------------------------|-------------------|
      | Peptide sequence | Sequence, Unique Peptide ID | peptide, sequence, Annotated Sequence, Peptide Sequence |
      | Precursor protein ID | Protein | Leading razor protein, UniProt ID, protein, Proteins, Protein, prot_acc, Accession, Master Protein Accessions, Protein ID |
      | Peptide start | start | start position, start, Start, pep_res_before, Positions in Master Proteins |
      | Peptide end | end | end position, end, End, pep_res_after, Positions in Master Proteins |
      | Samples intensity | user selected | user selected |
      | Modifications | Unique Peptide ID | Modified Sequence, Modifications, modified_peptide |

- **Functional annotations (optional):** MBPDB results or custom bioactivity databases
  - **Note:** MBPDB database searches must be performed using the web version at https://mbpdb.nws.oregonstate.edu/peptiline/
  - Download results as TSV and upload to local installation, or use example file: `examples/example_MBPDB_search.tsv`

    - Table 2. Required Columns for Functional Database Integration

      | Required columns | Explanation |
      |------------------|-------------|
      | search_peptide | matches your database sequences exactly—it's the key for linking functional and peptidomic data. |
      | peptide | Shows the database sequence found, which may differ slightly due to homology settings and partial matching. |
      | function | Contains functional data or other relevant information used for grouping the data. |

- **Protein sequences (optional):** FASTA files for custom protein variant analysis

See the `examples/` directory for sample datasets.


## Reproducibility

### Manuscript Case Study - Bitterness in Aged Cheddar Cheese
The manuscript case study analyzed bitter peptide profiles in aged cheddar cheese. Example files in the `examples/` directory provide representative data for testing the workflow orginating from the https://github.com/Kuhfeldrf/Bitterness-in-Aged-Cheddar-Cheese.

**Module 1: Data Transformation** (1-2 min)
- Upload: `example_peptide_data.csv`, `example_group_definition.json`, `example_MBPDB_search.tsv`, `example_fasta.fasta`
- Resolve multi-protein peptides
- Output: Merged dataset CSV

**Module 2: Descriptive Analysis** (30-60 sec)
- Upload merged dataset from Module 1
- Select parameters, generate visualizations

**Module 3: Heatmap Visualization** (20-40 sec per protein)
- Upload merged dataset, select protein, generate heatmap

**Total Runtime:** ~10-15 minutes


**Data Processing Workflow:**

**Module 1 (Data Transformation):**

   ***Step 1: Upload Data***
   - Peptidomic File: `examples/example_peptide_data.csv`
   - Functional Data File: `examples/example_MBPDB_search.tsv`
   - FASTA File: `examples/example_fasta.fasta`

   ***Step 2 (Optional): Organize Peptides with Multiple Protein Mappings***
   - Selected 'Yes' to process peptides mapped to multiple proteins
   - Applied the following decisions:

   | Combination | Occurrences | Protein ID | Decision |
   |-------------|-------------|------------|----------|
   | 1 | 1 | G9G9X6 | Remove |
   | | | P00711 | New |
   | 2 | 2 | A5D974 | Remove |
   | | | P02663 | New |
   | 3 | 4 | C6KGD7 | Remove |
   | | | C6KGD8 | New |
   | | | C6KGD9 | Remove |
   | 4 | 1 | A5D7J7 | Remove |
   | | | P02663 | New |
   | 5 | 20 | P31096 | New |
   | | | Q58DM6 | Remove |
   | 6 | 1 | A5D980 | Remove |
   | | | P02662 | New |
   | 7 | 6 | P02754 | New |
   | | | Q9BDG3 | Remove |
   | 8 | 928 | P02666A1 | Custom: P02666 |
   | | | P02666A2 | Remove |
   | 9 | 1 | P02754 | New |
   | | | Q9TRB9 | Remove |

   ***Step 3 (Optional): Assign Study Variables for Data Grouping***
   - Group Definitions File: `examples/example_group_definition.json`

   ***Step 4: Process & Export Data***
   - Clicked the Generate/Update Data button and exported the merged dataset as `examples/example_merged_dataframe_A1_and_A2_Beta_Casein.csv`
   - Additional processed data outputs were generated and are available in the `supplementals/` folder as Supplemental Tables 2-10

   ***Manual Edits to Aggregate β-casein Variants:***

   The file `examples/example_merged_dataframe_A1_and_A2_Beta_Casein.csv` was manually edited to aggregate β-casein variants and clean up protein names with the following changes:

   - **'Protein' and 'Positions in Proteins' columns:**
   - P02666A1 → P02666
   - P02666A2 → P02666
   - **'protein_name' column:**
   - CASB_BOVIN BetaA1-casein → CASB_BOVIN Beta-casein
   - CASB_BOVIN BetaA2-casein → CASB_BOVIN Beta-casein
   - Removed CAS*_BOVIN from all rows
   - These edits were saved as `examples/example_merged_dataframe.csv`

**Module 2 (Descriptive Analysis):**
   - Uploaded transformed dataset: `examples/example_merged_dataframe.csv`

   ***Supplemental Figures 1-5 Settings:***

   | Setting | Supp Fig 1 (SPLOM) | Supp Fig 2 (Functional Distribution) | Supp Fig 3 (Protein Origins) | Supp Fig 4 (Peptide Count) | Supp Fig 5 (Summed Absorbance) |
   |---------|-------------------|--------------------------------------|------------------------------|---------------------------|-------------------------------|
   | **Select Groups** | Threshold, Low, Moderate, Extreme | Non_bitter, Bitter | Non_bitter, Bitter | All (Threshold, Low, Moderate, Extreme, Non_bitter, Bitter) | All (Threshold, Low, Moderate, Extreme, Non_bitter, Bitter) |
   | **Select Proteins** | NA | NA | Beta-casein, AlphaS1-casein, AlphaS2-casein, Kappa-casein | NA | NA |
   | **Select Functions** | NA | All Functional Peptides | NA | NA | NA |
   | **Plot Filter** | No Filter | Selected Function(s) | Selected Protein(s) | No Filter | No Filter |
   | **Plot Type** | Corr Scatter Plots | Stacked Bar Plots | Stacked Bar Plots | Grouped Bar Plots | Grouped Bar Plots |
   | **Data Type** | Absorbance | Absorbance | Absorbance | Count | Absorbance |
   | **Scale Absorbance** | - | Absolute | Absolute | Absolute | Absolute |
   | **Plot Orientation** | - | By Sample | By Sample | By Sample | By Sample |
   | **Log10 Transform** | ✓ | - | - | - | - |
   | **Group Unselected** | - | - | ✓ | - | - |
   | **Color Scheme** | HSV | HSV | HSV | HSV | HSV |

**Module 3 (Heatmap Visualization):**
   - Upload Protein FASTA Files:
      - FASTA File: `examples/example_fasta.fasta`
      - ☐ Query UniProt for missing protein information 

   ***Figures 6A, 7A, and 7B Settings:***

   | Setting | Figure 6A (Bitter Samples) | Figure 7A (A1+A2 Merged) | Figure 7B (A1 & A2 Separated) |
   |---------|---------------------------|-------------------------|-------------------------------|
   | **Merged Data File** | `example_merged_dataframe.csv` | `example_merged_dataframe.csv` | `example_merged_dataframe_A1_and_A2_Beta_Casein.csv` |
   | **Update Labels** | - | - | βA1 Bitter, βA2 Bitter, βA1 Non-bitter, βA2 Non-bitter |
   | **Re-order Samples** | - | - | βA1 Bitter, βA2 Bitter, βA1 Non-bitter, βA2 Non-bitter |
   | **Select Groups** | Bitter | Bitter, Non_bitter | Bitter, Non_bitter |
   | **Select Proteins** | P02666 - CASB_BOVIN Beta-casein | P02666 - CASB_BOVIN Beta-casein | P02666A1 - BetaA1-casein, P02666A2 - BetaA2-casein |
   | **Plot Filter** | All Peptides | All Peptides | All Peptides |
   | **Plot Averaged Data** | yes | yes | yes |
   | **Plot Orientation** | Portrait | Landscape | Landscape |


## External Dependencies

### MBPDB (Milk Bioactive Peptide Database)
- **Access:** https://mbpdb.nws.oregonstate.edu (active, maintained)
- **⚠️ Local Limitation:** Database search NOT included in repository
  - Web version only: https://mbpdb.nws.oregonstate.edu/peptiline/
  - Local use: Upload pre-downloaded TSV files (see `examples/example_MBPDB_search.tsv`)
- **Reproducibility:** Manuscript version archived in `examples/` directory

### UniProt
- **Access:** https://www.uniprot.org (active, regularly updated)
- **Offline use:** Upload FASTA files directly (see `examples/example_fasta.fasta`)

### Python Packages
All dependencies in `notebook_requirements.txt` (standard PyPI packages). See [Dependencies](#dependencies).

## Documentation
- **Installation instructions:** This README
- **Module-specific guides:** See individual notebook interactive headers
- **Dependencies:** `notebook_requirements.txt`
- **Example data:** `examples/` directory
- **API utilities:** `utils/` directory

## Dependencies
Core dependencies include:
- Jupyter Notebook (≥7.4.5)
- pandas
- numpy
- matplotlib
- seaborn
- plotly
- scikit-learn
See `notebook_requirements.txt` for complete list with version specifications. Note: ipywidgets has been removed; notebooks with interactive widgets are deprecated.

## Project Structure

> This is the **legacy notebook layout** (v1.0). The production Django app is at https://github.com/mbpdb/mbpdb.

```
peptiline/  (legacy notebooks — v1.0)
├── data_transformation.ipynb      # Module 1: Data processing
├── data_analysis.ipynb            # Module 2: Descriptive analysis
├── heatmap_visualization.ipynb    # Module 3: Sequence visualization
├── _settings.py                   # Configuration file
├── notebook_requirements.txt      # Python dependencies
├── utils/                         # UniProt API utilities
├── examples/                      # Case study data
├── README.md                      # This file
└── LICENSE                        # MIT License
```

## Support and Contact
- **Bug reports and feature requests:** Open an issue on GitHub
- **General inquiries:** contact-mbpdb@oregonstate.edu
- **MBPDB Platform:** https://mbpdb.nws.oregonstate.edu

## Troubleshooting


### Timeouts or script errors in web browser
If you encounter timeouts or script errors while using PeptiLine in your web browser, try refreshing the page. This often resolves timeout issues or script execution problems.

### Memory issues with large datasets
If you encounter memory errors when processing large datasets:
- Close other applications to free up RAM
- Process data in smaller batches if possible
- Reduce the number of samples analyzed simultaneously
- Consider using a machine with more RAM for very large datasets (>10,000 peptides)
- Monitor memory usage using system tools during processing

### Connection timeouts 
If API calls to UniProt time out:
- Check your internet connection
- Use local files instead:
  - FASTA files for protein sequences
  - TSV files for functional annotations
- Retry the operation after a brief wait
- Note: Local installation does not allow direct MBPDB database searches. Use the web version at https://mbpdb.nws.oregonstate.edu/peptiline/ for MBPDB functionality, or upload pre-downloaded MBPDB results as TSV files

### Python version compatibility
If you experience import errors or unexpected behavior:
- Verify Python version: `python --version` (should be 3.12 or higher)
- Ensure virtual environment is activated
- Reinstall dependencies: `pip install -r notebook_requirements.txt --upgrade`

## Citation
**Note:** Manuscript under review at PLOS Computational Biology. Citation will be updated upon acceptance.

```
Kuhfeld R, Nielsen SD-H, Dallas DC. PeptiLine: An Interactive Platform for 
Customizable Functional Peptidomic Analysis. Submitted to PLOS Computational Biology. 2025.
```

**BibTeX:**
```bibtex
@article{kuhfeld2025peptiline,
  title={PeptiLine: An Interactive Platform for Customizable Functional Peptidomic Analysis},
  author={Kuhfeld, Russell and Nielsen, S{\o}ren D-H and Dallas, David C},
  journal={PLOS Computational Biology},
  year={2025},
  note={Submitted, under review}
}
```

## License
This project is licensed under the MIT License - an Open Source Initiative (OSI) compliant license. See the [LICENSE](LICENSE) file for details.

## Authors
- Russell Kuhfeld (Oregon State University)
- Søren D-H. Nielsen (Arla Foods Ingredients Group P/S)
- David C. Dallas (Oregon State University)


## Version History
- **v1.0** (2025) - Initial release accompanying PLOS Computational Biology submission (this repository)
- **v2.0** (2026) - Converted to Django web application; source at https://github.com/mbpdb/mbpdb

## Contributing
We welcome community contributions! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with a clear description of changes

For major changes, please open an issue first to discuss proposed modifications.

## Important Note: MBPDB Search Functionality

**⚠️ MBPDB database search is NOT included in local installations.**

**Why:** MBPDB database and search backend are proprietary components of the parent MBPDB application.

**Options:**
1. **Web version:** https://mbpdb.nws.oregonstate.edu/peptiline/ (full MBPDB search)
2. **Local:** Upload pre-downloaded MBPDB TSV files (see `examples/example_MBPDB_search.tsv`)
3. **Custom:** Upload your own functional annotations (TSV format, see Table 2)

**Available locally:** Data transformation, UniProt retrieval, all visualizations, MBPDB TSV uploads

---

## Software Archive
Manuscript version archived with journal submission (supplementary material) for long-term reproducibility.

**Archived:** Complete source code, example datasets, dependencies (`notebook_requirements.txt`), documentation, supplemental figures

**Version Control:** GitHub repository (https://github.com/Kuhfeldrf/peptiline), manuscript version tagged as v1.0

**Upon publication:** Will be updated with permanent archive reference (e.g., Zenodo DOI)
