"""
Group variable assignment logic.
Extracted from notebook GroupProcessing class (Cell 2).
"""
import json
import re
from collections import defaultdict

import pandas as pd

from .data_loader import get_filtered_columns


# ---------------------------------------------------------------------------
# Technical duplicate detection and collapsing
# ---------------------------------------------------------------------------

def _get_base_name(col):
    """Strip technical duplicate markers from a column name to get the biological replicate name."""
    # Suffix ( repN) parenthetical patterns: ( rep1), (rep2), ( rep 1), ( Rep2), etc.
    cleaned = re.sub(r'\s*\(\s*rep\s*\d+\s*\)$', '', col, flags=re.IGNORECASE)
    if cleaned != col:
        return cleaned.strip()

    # Suffix _repN underscore patterns: _rep1, _rep2, _rep_1, _Rep2, etc.
    cleaned = re.sub(r'[\s_\-]+rep[\s_\-]*\d+$', '', col, flags=re.IGNORECASE)
    if cleaned != col:
        return cleaned.strip()

    # Suffix _dup patterns: _dup, _dup1, _dup_1, _duplicate, _duplicates
    cleaned = re.sub(r'[\s_\-]+dup(?:licate[s]?)?[\s_\-]*\d*$', '', col, flags=re.IGNORECASE)
    if cleaned != col:
        return cleaned.strip()

    # Prefix rep patterns: rep1_, rep2_, rep_1_
    cleaned = re.sub(r'^rep[\s_\-]*\d+[\s_\-]+', '', col, flags=re.IGNORECASE)
    if cleaned != col:
        return cleaned.strip()

    # Prefix dup patterns: dup_, dup1_, duplicate_
    cleaned = re.sub(r'^dup(?:licate[s]?)?[\s_\-]+\d*', '', col, flags=re.IGNORECASE)
    if cleaned != col:
        return cleaned.strip()

    return col


def detect_technical_duplicates(columns):
    """
    Detect technical duplicates in a list of column names.

    Criteria:
    - A column name appears more than once (exact duplicate column names)
    - A column name ends with a ( repN) suffix: ( rep1), ( rep2), (rep1), etc.
    - A column name ends with a _repN suffix: _rep1, _rep2, _rep_1, etc.
    - A column name contains a _dup pattern: _dup, _dup1, _duplicate, etc.

    Returns:
        dict: {biological_replicate_name: [original_col_names]}
        Only includes groups where >1 columns map to the same base, or where
        the column name differs from its base (i.e. a rename is needed).
    """
    groups = defaultdict(list)
    for col in columns:
        groups[_get_base_name(col)].append(col)
    return {
        base: orig_cols
        for base, orig_cols in groups.items()
        if len(orig_cols) > 1 or orig_cols[0] != base
    }


def collapse_technical_duplicates(df, columns):
    """
    Average technical duplicate abundance columns into single biological replicate columns.

    Detects duplicates using detect_technical_duplicates() then collapses them
    by averaging. Uses positional indexing throughout to handle exact duplicate
    column names safely.

    Args:
        df: source DataFrame
        columns: list of abundance column names to check for duplicates

    Returns:
        (collapsed_df, dup_mapping, new_columns)
        - collapsed_df: DataFrame with technical duplicates averaged
        - dup_mapping: {bio_rep_name: [original_col_names]} (only dups/renames)
        - new_columns: ordered list of abundance column names after collapsing
    """
    dup_groups = detect_technical_duplicates(columns)
    if not dup_groups:
        return df.copy(), {}, list(columns)

    collapsed_df = df.copy()
    dup_mapping = {}

    for base, orig_cols in dup_groups.items():
        cur_cols = list(collapsed_df.columns)
        used_positions = set()
        series_list = []

        # Positional indexing: each occurrence in orig_cols gets a distinct position
        for col_name in orig_cols:
            for pos, c in enumerate(cur_cols):
                if c == col_name and pos not in used_positions:
                    series_list.append(pd.to_numeric(collapsed_df.iloc[:, pos], errors='coerce'))
                    used_positions.add(pos)
                    break

        if not series_list:
            continue

        avg = pd.concat(series_list, axis=1).mean(axis=1) if len(series_list) > 1 else series_list[0]

        # Rebuild df: drop all original positions, insert the averaged column at the first position
        keep = [i for i in range(len(cur_cols)) if i not in used_positions]
        insert_at = min(used_positions, default=len(keep))
        insert_at = min(insert_at, len(keep))
        temp = collapsed_df.iloc[:, keep].copy()
        temp.insert(insert_at, base, avg.values)
        collapsed_df = temp

        dup_mapping[base] = list(orig_cols)

    # Build ordered new column list (first occurrence wins)
    seen = set()
    new_columns = []
    for col in columns:
        base = _get_base_name(col)
        target = base if base in dup_mapping else col
        if target not in seen:
            new_columns.append(target)
            seen.add(target)

    return collapsed_df, dup_mapping, new_columns


def apply_tech_dup_mapping(df, tech_dup_mapping):
    """
    Apply a pre-computed technical duplicate mapping to a DataFrame.

    Used during final processing to collapse technical duplicates in the
    peptidomic data (pd_results / pd_results_cleaned) before group averaging.

    Args:
        df: source DataFrame
        tech_dup_mapping: {base_name: [original_col_names]}

    Returns:
        DataFrame with technical duplicates collapsed
    """
    if not tech_dup_mapping:
        return df.copy()

    collapsed_df = df.copy()
    for base, orig_cols in tech_dup_mapping.items():
        cur_cols = list(collapsed_df.columns)
        used_positions = set()
        series_list = []

        for col_name in orig_cols:
            for pos, c in enumerate(cur_cols):
                if c == col_name and pos not in used_positions:
                    series_list.append(pd.to_numeric(collapsed_df.iloc[:, pos], errors='coerce'))
                    used_positions.add(pos)
                    break

        if not series_list:
            continue

        avg = pd.concat(series_list, axis=1).mean(axis=1) if len(series_list) > 1 else series_list[0]

        keep = [i for i in range(len(cur_cols)) if i not in used_positions]
        insert_at = min(used_positions, default=len(keep))
        insert_at = min(insert_at, len(keep))
        temp = collapsed_df.iloc[:, keep].copy()
        temp.insert(insert_at, base, avg.values)
        collapsed_df = temp

    return collapsed_df


def compute_bio_rep_columns(raw_columns, tech_dup_mapping):
    """
    Compute the post-collapse column list given raw columns and a tech dup mapping.

    Each tech rep group's columns are replaced by a single bio rep name (inserted at
    the position of the first tech rep in the list).

    Args:
        raw_columns: list of column names before collapse
        tech_dup_mapping: {bio_rep_name: [tech_rep_col_names]}

    Returns:
        list of column names with tech rep groups collapsed to their bio rep names
    """
    col_to_bio_rep = {}
    for bio_rep, tech_rep_cols in tech_dup_mapping.items():
        for tc in tech_rep_cols:
            col_to_bio_rep[tc] = bio_rep

    result = []
    bio_rep_inserted = set()
    for col in raw_columns:
        bio_rep = col_to_bio_rep.get(col)
        if bio_rep is not None:
            if bio_rep not in bio_rep_inserted:
                result.append(bio_rep)
                bio_rep_inserted.add(bio_rep)
            # else: skip — already inserted as part of the same bio rep group
        else:
            result.append(col)
    return result


def parse_group_json(json_content, available_columns):
    """
    Parse a group definition JSON file.

    Args:
        json_content: string or bytes of JSON content
        available_columns: list of column names available in the dataset

    Returns:
        tuple: (group_data dict, error_message or None)
            group_data format: {
                "1": {"grouping_variable": "name", "abundance_columns": ["col1", "col2"]},
                ...
            }
    """
    try:
        if isinstance(json_content, bytes):
            json_content = json_content.decode('utf-8')
        simplified_data = json.loads(json_content)
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        return None, f"Invalid JSON file: {str(e)}"

    # Check for missing columns
    available_set = set(available_columns)
    missing_columns = []

    for group_name, columns in simplified_data.items():
        if isinstance(columns, list):
            for column in columns:
                if column not in available_set:
                    missing_columns.append(column)

    if missing_columns:
        missing_unique = sorted(set(missing_columns))
        return None, f"The following columns are not present in the current dataset: {', '.join(missing_unique)}"

    # Use the group name directly as the key (no numeric wrapping)
    group_data = {}
    for group_name, abundance_cols in simplified_data.items():
        if isinstance(abundance_cols, list):
            group_data[group_name] = {
                'grouping_variable': group_name,
                'abundance_columns': abundance_cols
            }
        else:
            return None, f"Invalid format: group '{group_name}' should map to a list of column names"

    return group_data, None


def build_group_data(selected_columns, group_name, existing_group_data=None):
    """
    Add a new group to the group data.

    Args:
        selected_columns: list of column names for this group
        group_name: name for the grouping variable
        existing_group_data: existing group_data dict to add to

    Returns:
        tuple: (updated group_data, error_message or None)
    """
    if not group_name:
        return existing_group_data, "Please enter a group name"
    if not selected_columns:
        return existing_group_data, "Please select at least one column"

    if existing_group_data is None:
        existing_group_data = {}

    existing_group_data[group_name] = {
        'grouping_variable': group_name,
        'abundance_columns': list(selected_columns)
    }

    return existing_group_data, None


def build_no_group_data(selected_columns):
    """
    Create individual groups for each selected column (no grouping).

    Args:
        selected_columns: list of column names

    Returns:
        tuple: (group_data, error_message or None)
    """
    if not selected_columns:
        return None, "Please select at least one abundance column"

    group_data = {}
    for column in selected_columns:
        group_data[column] = {
            'grouping_variable': column,
            'abundance_columns': [column]
        }

    return group_data, None


def validate_group_data(group_data, df):
    """
    Validate that group data references valid columns in the dataframe.

    Args:
        group_data: group_data dict
        df: pandas DataFrame

    Returns:
        tuple: (is_valid, error_message or None)
    """
    if not group_data:
        return False, "No group data provided"

    df_columns = set(df.columns)
    for group_id, group_info in group_data.items():
        missing = [col for col in group_info['abundance_columns'] if col not in df_columns]
        if missing:
            return False, f"Group '{group_info['grouping_variable']}' references missing columns: {', '.join(missing)}"

    return True, None


def export_group_data_json(group_data):
    """
    Export group data to simplified JSON format for download.

    Returns:
        str: JSON string
    """
    simplified = {}
    for _, group_info in group_data.items():
        simplified[group_info['grouping_variable']] = group_info['abundance_columns']
    return json.dumps(simplified, indent=4)
