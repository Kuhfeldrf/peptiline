# Data Validation Utilities

This folder contains utilities for validating and processing data used across the MBPDB notebooks. The validation system uses Pandera for robust DataFrame validation and custom validators for other file types.

## Dependencies

```
pandas
pandera
```

## Files

- `data_validators.py`: Contains validator classes for different types of input data
- `data_schemas.py`: Defines Pandera schemas and other validation schemas
- `data_utils.py`: Provides utility functions for data loading and validation

## Usage

### Basic Usage

```python
from utils.data_utils import (
    load_and_validate_peptide_data,
    load_and_validate_fasta,
    load_and_validate_group_definition,
    load_and_validate_merged_data
)

# Load and validate peptide data (uses Pandera schema)
peptide_df = load_and_validate_peptide_data('path/to/peptide_data.csv')

# Load and validate FASTA file
sequences = load_and_validate_fasta('path/to/sequences.fasta')

# Load and validate group definitions
groups = load_and_validate_group_definition('path/to/groups.json')

# Load and validate merged data (uses Pandera schema)
merged_df = load_and_validate_merged_data('path/to/merged_data.csv')
```

### Advanced Usage with Pandera Schemas

You can use the Pandera schemas directly for more control:

```python
from utils.data_schemas import PeptideDataSchema, MergedDataSchema
import pandas as pd

# Validate DataFrame using schema
df = pd.read_csv('path/to/data.csv')
validated_df = PeptideDataSchema.validate(df)

# Use schema for type hints and runtime checking
@pa.check_types
def process_peptide_data(df: PeptideDataSchema) -> PeptideDataSchema:
    # Your processing logic here
    return df
```

### Direct Validator Usage

For more control over the validation process:

```python
from utils.data_validators import PeptideDataValidator

validator = PeptideDataValidator()
is_valid = validator.validate('path/to/peptide_data.csv')

if not is_valid:
    print("Validation errors:")
    for error in validator.get_errors():
        print(f"- {error}")
```

## Data Types and Requirements

### Peptide Data (CSV)
- Validated using Pandera schema
- Required columns: sequence, name
- Sequence validation: only valid amino acid characters
- No null values allowed

### FASTA Files
- Must start with '>' header line
- Valid amino acid sequences
- Header format validation
- Standard FASTA format

### Group Definition (JSON)
- Must be a valid JSON dictionary
- All group members must be lists of strings
- Required groups: Threshold, Low, Moderate, Extreme, Non_bitter, Bitter

### Merged Data (CSV)
- Validated using Pandera schema
- Required columns: sequence, name
- No null values allowed
- Additional validations based on specific requirements

## Error Handling

All validation functions return `None` and print error messages when validation fails. This allows for easy error handling in notebooks:

```python
df = load_and_validate_peptide_data('path/to/data.csv')
if df is None:
    print("Failed to load data, please check the errors above")
else:
    # Continue with processing
    pass
```