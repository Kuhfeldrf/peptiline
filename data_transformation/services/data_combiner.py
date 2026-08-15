"""
Merge peptidomic + functional data, calculate averages.
Extracted from notebook CombineAverageDataframes class (Cell 5).
"""
import warnings

import pandas as pd
import numpy as np

from .data_loader import find_species, extract_protein_id


def add_protein_info(df, protein_dict):
    """Add protein species and name columns based on Protein column."""
    df = df.copy()

    if 'protein_species' not in df.columns:
        df['protein_species'] = 'Unknown'
    if 'protein_name' not in df.columns:
        df['protein_name'] = 'Unknown'

    for idx, row in df.iterrows():
        if isinstance(row['Protein'], list):
            proteins = [str(p) for p in row['Protein']]
        else:
            protein_str = str(row['Protein'])
            if ';' in protein_str:
                proteins = protein_str.split(';')
            elif '/' in protein_str:
                proteins = protein_str.split('/')
            elif ',' in protein_str:
                proteins = protein_str.split(',')
            else:
                proteins = [protein_str]

        if proteins and proteins[0] != '' and proteins[0] != 'nan':
            protein_raw = proteins[0].strip()

            if '|' in protein_raw:
                parts = protein_raw.split('|')
                protein_id = parts[1].strip() if len(parts) >= 2 else protein_raw
            else:
                protein_id = protein_raw

            protein_info = protein_dict.get(protein_id, {})

            if not protein_info and '|' in protein_raw:
                parts = protein_raw.split('|')
                if len(parts) == 3:
                    name_species = parts[2]
                    if '_' in name_species:
                        name, species_abbrev = name_species.split('_', 1)
                        species = find_species(species_abbrev)
                        protein_info = {'name': name if name else protein_id, 'species': species}
                        protein_dict[protein_id] = protein_info

            if df.at[idx, 'protein_species'] == 'Unknown':
                df.at[idx, 'protein_species'] = protein_info.get('species', "Unknown")
            if df.at[idx, 'protein_name'] == 'Unknown':
                df.at[idx, 'protein_name'] = protein_info.get('name', "Unknown")

    # Reorder columns
    all_cols = list(df.columns)
    remaining_cols = [col for col in all_cols if col not in ['protein_species', 'protein_name']]
    if 'Protein' in remaining_cols:
        insert_pos = remaining_cols.index('Protein') + 1
    else:
        insert_pos = 0
    new_cols = remaining_cols[:insert_pos] + ['protein_species', 'protein_name'] + remaining_cols[insert_pos:]
    return df.reindex(columns=new_cols)


def extract_bioactive_peptides(mbpdb_results):
    """Extract the list of bioactive peptide matches from MBPDB search."""
    if mbpdb_results is None or mbpdb_results.empty:
        return None, None

    mbpdb_results_cleaned = mbpdb_results.copy()
    mbpdb_results_cleaned.dropna(subset=['search_peptide'], inplace=True)
    if 'protein_id' in mbpdb_results_cleaned.columns:
        mbpdb_results_cleaned = mbpdb_results_cleaned[mbpdb_results_cleaned['protein_id'] != 'None']

    available_columns = mbpdb_results_cleaned.columns.tolist()
    agg_dict = {}

    if 'peptide' in available_columns:
        agg_dict['peptide'] = 'first'
    if 'protein_id' in available_columns:
        agg_dict['protein_id'] = 'first'

    optional_columns = {
        'protein_description': 'first', '% Alignment': 'first',
        'species': 'first', 'intervals': 'first',
        'additional_details': 'first', 'ic50': 'first',
        'inhibition_type': 'first', 'inhibited_microorganisms': 'first',
        'ptm': 'first', 'title': 'first', 'authors': 'first',
        'abstract': 'first', 'doi': 'first',
        'search_type': 'first', 'scoring_matrix': 'first'
    }

    for col, agg_func in optional_columns.items():
        if col in available_columns:
            agg_dict[col] = agg_func

    if 'function' in available_columns:
        agg_dict['function'] = lambda x: list(x.dropna().unique())

    if not agg_dict:
        return mbpdb_results_cleaned, mbpdb_results_cleaned

    try:
        mbpdb_results_grouped = mbpdb_results_cleaned.groupby('search_peptide').agg(agg_dict).reset_index()
        if 'function' in mbpdb_results_grouped.columns:
            mbpdb_results_grouped['function'] = mbpdb_results_grouped['function'].apply(
                lambda x: '; '.join([str(func) for func in x if str(func) != 'nan']) if isinstance(x, list) else str(x)
            )
        return mbpdb_results_cleaned, mbpdb_results_grouped
    except Exception:
        return mbpdb_results_cleaned, mbpdb_results_cleaned


def create_unique_id(row):
    """Creates a Unique Peptide ID for each peptide row."""
    sequence = row['Sequence']
    if isinstance(sequence, list):
        sequence = ','.join(sequence)
    else:
        sequence = str(sequence).strip()

    for mod_col in ['Peptide Modified Sequence', 'Modified sequence', 'Modified Sequence',
                     'EG.ModifiedPeptide', 'modified_peptide']:
        if mod_col in row and pd.notna(row[mod_col]):
            unique_id = str(row[mod_col]).strip()
            break
    else:
        if 'Modifications' in row and pd.notna(row['Modifications']):
            unique_id = sequence + "_" + str(row['Modifications']).strip()
        else:
            unique_id = sequence

    return str(unique_id).strip().rstrip('_')


_PLACEHOLDER_PROTEIN_VALUES = {'Unknown', 'unknown', 'UNKNOWN', '', ' ', 'None', 'none', 'nan', 'NaN'}


def _row_signature(row):
    """NaN-safe, order-insensitive-for-lists string signature of a row's
    values, used to test whether rows are true duplicates (not just sharing
    a peptide identity)."""
    parts = []
    for v in row:
        if isinstance(v, (list, tuple, np.ndarray)):
            parts.append('|'.join(sorted(str(x) for x in v)))
        elif v is None or (pd.api.types.is_scalar(v) and pd.isna(v)):
            parts.append('\x00')
        else:
            parts.append(str(v))
    return '\x1f'.join(parts)


def _collect_proteins(values):
    """Join Protein/Positions-in-Proteins values from multiple rows (each
    possibly a list, or a ';'-separated string) into one deduplicated,
    order-preserving '; '-joined string, matching how Proteome Discoverer,
    MaxQuant, PEAKS, and Spectronaut already represent multi-protein
    peptides in a single row."""
    seen = []
    for v in values:
        if v is None or (not isinstance(v, (list, tuple, np.ndarray)) and pd.isna(v)):
            continue
        items = v if isinstance(v, (list, tuple, np.ndarray)) else str(v).split(';')
        for item in items:
            item = str(item).strip()
            if item and item not in _PLACEHOLDER_PROTEIN_VALUES and item not in seen:
                seen.append(item)
    return '; '.join(seen) if seen else 'Unknown'


def _collapse_multiprotein_duplicate_rows(df):
    """Collapse rows that share a 'Unique Peptide ID' but differ only in
    'Protein' / 'Positions in Proteins' into a single row with those columns
    joined — the layout Proteome Discoverer, MaxQuant, PEAKS, and
    Spectronaut already use for a peptide mapped to multiple proteins.

    Skyline instead reports one row per (peptide, protein) pair with
    otherwise identical data (abundances included), so without this the
    same peptide's abundance would silently persist as several rows sharing
    one 'Unique Peptide ID' — inflating any sum-based downstream analysis
    and contradicting the "avoids errors from duplicate peptides" guarantee
    the ID is meant to provide (Table 1, footnote a).

    Only collapses when every non-protein column is identical (NaN-safe)
    across the group; groups that differ elsewhere are left untouched since
    that indicates genuinely distinct measurements, not a protein-mapping
    artifact, and must not be silently merged.
    """
    if 'Unique Peptide ID' not in df.columns:
        return df

    uid = df['Unique Peptide ID']
    dup_mask = uid.notna() & uid.duplicated(keep=False)
    if not dup_mask.any():
        return df

    # Per-protein descriptor columns some engines report alongside 'Protein'
    # (e.g. Skyline's 'Protein Name'/'Protein Description'/'Protein
    # Accession', which duplicate 'Protein' one-for-one) — these legitimately
    # differ across a peptide's protein rows and must be joined rather than
    # required to match for a group to be considered collapsible.
    protein_like_cols = [c for c in (
        'Protein', 'Positions in Proteins', 'Protein Name',
        'Protein Description', 'Protein Accession',
        'Master Protein Descriptions', 'Leading razor protein',
    ) if c in df.columns]
    other_cols = [c for c in df.columns if c not in protein_like_cols + ['Unique Peptide ID']]

    dup_df = df.loc[dup_mask]
    signatures = dup_df[other_cols].apply(_row_signature, axis=1) if other_cols else pd.Series('', index=dup_df.index)

    collapsible_ids = [
        pid for pid, sig_group in signatures.groupby(dup_df['Unique Peptide ID'])
        if sig_group.nunique() == 1
    ]
    if not collapsible_ids:
        return df

    collapsible_mask = dup_mask & uid.isin(collapsible_ids)
    keep_asis = df.loc[~collapsible_mask]
    to_collapse = df.loc[collapsible_mask]

    agg = {c: 'first' for c in other_cols}
    agg.update({c: _collect_proteins for c in protein_like_cols})

    collapsed = to_collapse.groupby('Unique Peptide ID', as_index=False).agg(agg)

    result = pd.concat([keep_asis, collapsed], ignore_index=True, sort=False)
    return result[df.columns]


def process_pd_results(pd_results_cleaned, mbpdb_results_grouped, protein_dict):
    """Process peptidomic results and merge with MBPDB data."""
    df = pd_results_cleaned.copy()

    if 'Protein' not in df.columns:
        df['Protein'] = pd.NA
    if 'Positions in Proteins' not in df.columns:
        df['Positions in Proteins'] = pd.NA

    df['Protein'] = df['Protein'].fillna('Unknown')
    df['Positions in Proteins'] = df['Positions in Proteins'].fillna('Unknown')

    if 'Sequence' not in df.columns:
        df['Sequence'] = pd.NA
        if 'Annotated Sequence' in df.columns:
            def extract_seq(annotated_seq):
                if pd.isna(annotated_seq):
                    return pd.NA
                if '.' in str(annotated_seq):
                    parts = str(annotated_seq).split('.')
                    if len(parts) > 1:
                        return parts[1]
                return annotated_seq
            df['Sequence'] = df['Annotated Sequence'].apply(extract_seq)

    df['Unique Peptide ID'] = df.apply(create_unique_id, axis=1)
    df = _collapse_multiprotein_duplicate_rows(df)

    # Extract start/stop positions
    try:
        for col in ['start', 'end']:
            if col not in df.columns:
                df[col] = pd.NA

        valid_position_mask = (
            ~df['Positions in Proteins'].str.contains(';', na=False) &
            ~df['Positions in Proteins'].str.contains('/', na=False) &
            (df['Positions in Proteins'] != 'Unknown')
        )
        single_positions = df.loc[valid_position_mask, 'Positions in Proteins']
        if not single_positions.empty:
            extracted = single_positions.str.extract(r'\[(\d+)-(\d+)\]')
            df.loc[valid_position_mask, 'start'] = pd.to_numeric(extracted[0], errors='coerce')
            df.loc[valid_position_mask, 'end'] = pd.to_numeric(extracted[1], errors='coerce')

        df['start'] = df['start'].astype('Int64')
        df['end'] = df['end'].astype('Int64')
    except Exception:
        pass

    # Reorder columns
    remaining_cols = [col for col in df.columns
                     if col not in ['Unique Peptide ID', 'Sequence', 'Protein',
                                  'Positions in Proteins', 'start', 'end']]
    columns_order = ['Unique Peptide ID', 'Sequence', 'Protein',
                    'Positions in Proteins', 'start', 'end'] + remaining_cols
    df = df[columns_order]

    # Merge with MBPDB results
    if mbpdb_results_grouped is not None and not mbpdb_results_grouped.empty:
        # Drop any pre-existing function column to avoid pandas merge suffix conflicts;
        # the mbpdb function column takes priority.
        df_to_merge = df.drop(columns=['function'], errors='ignore')
        merged_df = pd.merge(df_to_merge, mbpdb_results_grouped,
                            right_on='search_peptide', left_on='Unique Peptide ID', how='left')

        # Handle comma-separated IDs
        comma_mask = merged_df['Unique Peptide ID'].str.contains(',', na=False)
        comma_rows = merged_df[comma_mask].copy()
        for idx, row in comma_rows.iterrows():
            unique_ids = row['Unique Peptide ID'].split(',')
            matches = mbpdb_results_grouped[mbpdb_results_grouped['search_peptide'].isin(unique_ids)]
            if not matches.empty:
                match = matches.iloc[0]
                for col in mbpdb_results_grouped.columns:
                    merged_df.loc[idx, col] = match[col]
    else:
        merged_df = df.copy()
        if 'function' not in df.columns:
            # No MBPDB results and no embedded function column — add empty placeholder
            merged_df['function'] = np.nan
        # else: preserve the existing function column from the uploaded file

    final_column_order = columns_order + [col for col in merged_df.columns if col not in columns_order]
    merged_df = merged_df[final_column_order]

    return merged_df


def calculate_group_abundance_averages(df, group_data):
    """Calculate group abundance averages."""
    # Check if all average columns already exist
    all_exist = all(
        f"Avg_{details['grouping_variable']}" in df.columns
        for details in group_data.values()
    )
    if all_exist:
        return df

    average_columns = {}
    for group_number, details in group_data.items():
        grouping_variable = details['grouping_variable']
        abundance_columns = details['abundance_columns']

        for col in abundance_columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        average_column_name = f"Avg_{grouping_variable}"
        average_columns[average_column_name] = df[abundance_columns].mean(axis=1, skipna=True)

    df = pd.concat([df, pd.DataFrame(average_columns)], axis=1)
    return df


def update_column_names_with_groups(df, group_data):
    """Update column names by adding grouping information."""
    column_to_groups = {}
    for group_number, details in group_data.items():
        grouping_variable = details['grouping_variable']
        for column in details['abundance_columns']:
            if column not in column_to_groups:
                column_to_groups[column] = []
            column_to_groups[column].append(grouping_variable)

    df_renamed = df.copy()
    rename_dict = {}
    for column in column_to_groups:
        if column in df.columns:
            groups_str = "; ".join(column_to_groups[column])
            new_name = f"{column} 'Grouped: ({groups_str})'"
            rename_dict[column] = new_name

    df_renamed = df_renamed.rename(columns=rename_dict)
    return df_renamed


def process_data(pd_results, pd_results_cleaned, mbpdb_results, group_data, protein_dict):
    """
    Main data processing function.

    Args:
        pd_results: original peptidomic DataFrame
        pd_results_cleaned: cleaned peptidomic DataFrame (after protein mapping)
        mbpdb_results: MBPDB search results DataFrame
        group_data: group definitions dict
        protein_dict: protein information dict

    Returns:
        pd.DataFrame or None
    """
    if pd_results is None or pd_results.empty:
        return None

    step = 'initializing'
    try:
        step = 'extracting bioactive peptides'
        mbpdb_results_cleaned, mbpdb_results_grouped = extract_bioactive_peptides(mbpdb_results)

        if pd_results_cleaned is None or pd_results_cleaned.empty:
            pd_results_cleaned = pd_results.copy()

        step = 'merging peptidomic results with MBPDB data'
        merged_df = process_pd_results(pd_results_cleaned, mbpdb_results_grouped, protein_dict)

        if group_data:
            step = 'calculating group abundance averages'
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                merged_df = calculate_group_abundance_averages(merged_df, group_data)
                step = 'renaming columns with group labels'
                merged_df = update_column_names_with_groups(merged_df, group_data)

        step = 'adding protein species/name info'
        final_df = add_protein_info(merged_df, protein_dict)

        step = 'cleaning protein IDs'
        if 'Protein' in final_df.columns:
            final_df['Protein'] = final_df['Protein'].apply(
                lambda x: extract_protein_id(x, protein_dict)
            )
            # extract_protein_id returns a list for multi-protein entries;
            # normalize to '; '-separated strings for consistent downstream use.
            final_df['Protein'] = final_df['Protein'].apply(
                lambda x: '; '.join(str(p) for p in x) if isinstance(x, list) else x
            )

        step = 'cleaning up placeholder values'
        for col in ['Protein', 'Positions in Proteins', 'protein_name', 'protein_species']:
            if col in final_df.columns:
                final_df[col] = final_df[col].replace(
                    ['Unknown', 'unknown', 'UNKNOWN', '', ' ', None, 'None', 'none', 'nan', 'NaN'],
                    pd.NA
                )

        return final_df
    except Exception as exc:
        raise RuntimeError(
            f'Data processing failed at step "{step}": {type(exc).__name__}: {exc}'
        ) from exc
