# UniProt Integration for Peptide Visualization Notebooks

This module provides integration with the UniProt database for the peptide visualization notebooks. It allows you to:

1. Automatically fetch missing protein information from UniProt
2. Search for and import specific proteins by UniProt ID
3. Import proteins directly from UniProt instead of uploading FASTA files

## Available Notebooks

We provide two notebooks that demonstrate the UniProt integration:

1. `data_transformation.ipynb` - Integrates UniProt with the data transformation notebook
2. `heatmap_visualization.ipynb` - Integrates UniProt with the heatmap visualization notebook

## How to Use

### Method 1: Use the Pre-integrated Notebooks

Simply run one of the integrated notebooks mentioned above. They include all the necessary code to enable UniProt integration.

### Method 2: Add Integration to Existing Notebooks

Add these lines to your existing notebook:

```python
# Import the UniProt integration
from utils.uniprot_integration import UniProtIntegration

# Create the integration with your data transformer
uniprot_integration = UniProtIntegration(data_transformer)

# Add the UniProt UI to the notebook
uniprot_integration.add_uniprot_ui()
```

## Features

### Automatic Detection of Missing Proteins

When you upload a merged data file, the integration will automatically check for proteins that are referenced in your data but not available in your uploaded FASTA files. It will then provide options to fetch these proteins from UniProt.

### Search for Specific Proteins

You can search for specific proteins by their UniProt ID (e.g., P01308 for human insulin) and import them directly into your visualization.

### Batch Import of Proteins

You can import multiple proteins at once either by:
1. Using the "Fetch All Missing" button to import all missing proteins detected in your data
2. Programmatically importing a list of UniProt IDs using the provided functions

## API Reference

### UniProtIntegration Class

The main class for integrating UniProt with the peptide visualization notebooks.

#### Methods:

- `__init__(data_transformer)` - Initialize with a reference to your data transformer
- `add_uniprot_ui([output_area])` - Add the UniProt UI components to the notebook. Optionally accepts a custom output area.
- `check_missing_proteins()` - Check for proteins in the merged data that are missing from the protein dictionary
- `fetch_batch_from_uniprot()` - Fetch all missing proteins from UniProt in batch mode

### Utility Functions

- `add_uniprot_integration_to_dataframe(df, uniprot_client=None)` - Add protein name and species information to a DataFrame by querying UniProt. The DataFrame must have a 'Master Protein Accessions' column.

## Requirements

- The `uniprot_client.py` module (included in utils)
- Internet connection to access the UniProt database
- Python packages: requests, pandas, ipywidgets

## Example Usage

```python
# Import the module
from utils.uniprot_integration import UniProtIntegration

# Initialize with a reference to your data transformer
uniprot_integration = UniProtIntegration(data_transformer)

# Add the UI components (optionally with a custom output area)
uniprot_integration.add_uniprot_ui()  # or add_uniprot_ui(output_area=custom_output)

# The integration will automatically check for missing proteins when enabled
# You can also manually check:
uniprot_integration.check_missing_proteins()

# To fetch all missing proteins at once (after enabling UniProt integration):
uniprot_integration.fetch_batch_from_uniprot()
```

## Troubleshooting

### Common Issues:

1. **Cannot connect to UniProt**: Make sure you have a working internet connection. UniProt servers may also occasionally be down for maintenance.

2. **UniProt ID not found**: Check that you are using the correct UniProt accession ID (e.g., P01308).

3. **No sequence returned**: Some proteins in UniProt may not have sequence information available.

4. **Slow performance**: Batch requests can be slow, especially for large numbers of proteins. Consider fetching only the proteins you need. 