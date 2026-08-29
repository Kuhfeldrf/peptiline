"""
Export logic: generate Excel/CSV exports, correlations.
Extracted from notebook ExportManager class (Cell 6).
"""
import io
import json
from collections import defaultdict
from itertools import combinations

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr


def _apply_scientific_notation(worksheet):
    """Format cells with |value| > 10,000 as scientific notation in an openpyxl worksheet."""
    for row in worksheet.iter_rows(min_row=2):  # Skip header row
        for cell in row:
            if isinstance(cell.value, (int, float)) and cell.value is not None:
                try:
                    if abs(cell.value) > 10000:
                        cell.number_format = '0.00E+00'
                except (TypeError, ValueError):
                    pass


def calculate_correlation(x, y, correlation_type='Pearson'):
    """Calculate correlation and return (correlation, p_value_string)."""
    if correlation_type == 'Pearson':
        corr, p = pearsonr(x, y)
    else:
        corr, p = spearmanr(x, y)

    if p < 0.001:
        p_str = "p < 0.001"
    else:
        p_str = f"p = {p:.3g}"

    return corr, p_str


def prepare_data(data, log_transform=True, log_base=10):
    """Prepare data based on log transform setting.

    ``log_base`` selects the logarithm base when ``log_transform`` is truthy
    (10 or 2). Log10 and Log2 are mutually exclusive in the UI.
    """
    if log_transform:
        return np.log2(data) if int(log_base) == 2 else np.log10(data)
    return data


def log_label(log_transform, log_base=10):
    """Human-readable name of the active transform, for export annotations."""
    if not log_transform:
        return 'None'
    return 'Log2' if int(log_base) == 2 else 'Log10'


def _write_summary_with_transform(writer, summary_df, transform_label,
                                  correlation_type, sheet_name='Summary'):
    """Write a correlation Summary sheet, prefaced with the parameters used.

    Records the correlation method and the log transform (Log10 / Log2 / None)
    that produced the correlations so the export is self-documenting.
    """
    meta = pd.DataFrame({
        'Parameter': ['Correlation Method', 'Log Transform'],
        'Value': [correlation_type, transform_label],
    })
    meta.to_excel(writer, sheet_name=sheet_name, index=False, startrow=0)
    summary_df.to_excel(writer, sheet_name=sheet_name, index=False,
                        startrow=len(meta) + 2)


def _convert_group_data_dict(merged_df, group_data):
    """
    Strip 'Grouped:' suffixes from column names and update group_data references.

    Returns:
        tuple: (updated_group_data, df_with_renamed_columns)
    """
    df = merged_df.copy()

    grouped_columns = [col for col in df.columns if " 'Grouped:" in str(col)]
    renamed_columns = {}
    for col in grouped_columns:
        base_col_name = col.split(" 'Grouped:")[0].strip()
        renamed_columns[col] = base_col_name
    df = df.rename(columns=renamed_columns)

    updated_group_data = {}
    for key, value in group_data.items():
        grouping_variable = value['grouping_variable']
        abundance_columns = value['abundance_columns']
        updated_abundance_columns = [col for col in abundance_columns if col in df.columns]
        updated_group_data[key] = {
            'grouping_variable': grouping_variable,
            'abundance_columns': updated_abundance_columns
        }

    return updated_group_data, df


def _get_abundance_cols(group_data):
    """Get the average abundance column names from group_data."""
    abundance_cols = []
    selected_groups = []
    for _, value in group_data.items():
        grouping_variable = value['grouping_variable']
        abundance_columns = value['abundance_columns']
        selected_groups.append(grouping_variable)
        if grouping_variable != abundance_columns:
            abundance_cols.append(f'Avg_{grouping_variable}')
        else:
            abundance_cols.append(abundance_columns)
    return abundance_cols, selected_groups


# ---------------------------------------------------------------------------
# Export: functional annotation results (TSV)
# ---------------------------------------------------------------------------

def export_mbpdb_results(mbpdb_df):
    """Export functional annotation results as TSV bytes."""
    if mbpdb_df is None or 'function' not in mbpdb_df.columns:
        return None
    return mbpdb_df.to_csv(sep='\t', index=False).encode('utf-8')


# ---------------------------------------------------------------------------
# Export: Group definitions (JSON)
# ---------------------------------------------------------------------------

def export_group_definitions(group_data):
    """Export group data as JSON bytes."""
    if not group_data:
        return None
    simplified = {}
    for _, group_info in group_data.items():
        simplified[group_info['grouping_variable']] = group_info['abundance_columns']
    return json.dumps(simplified, indent=4).encode('utf-8')


# ---------------------------------------------------------------------------
# Export: Merged dataset (CSV)
# ---------------------------------------------------------------------------

def export_merged_dataset(merged_df):
    """Export merged dataframe as CSV bytes."""
    if merged_df is None:
        return None
    return merged_df.to_csv(index=False).encode('utf-8')


# ---------------------------------------------------------------------------
# Export: Sequence list (CSV)
# ---------------------------------------------------------------------------

def extract_sequences_by_name(merged_df, group_data):
    """Extract unique peptide IDs per group based on non-zero average abundance."""
    result_dict = {}
    available_columns = merged_df.columns.tolist()

    for group_id, group_info in group_data.items():
        grouping_variable = f"Avg_{group_info['grouping_variable']}"
        if grouping_variable in available_columns:
            mask = merged_df[grouping_variable] > 0
            unique_peptides = merged_df.loc[mask, 'Unique Peptide ID'].unique().tolist()
            result_dict[group_info['grouping_variable']] = unique_peptides

    if not result_dict:
        return pd.DataFrame()

    max_length = max(len(v) for v in result_dict.values())
    padded_dict = {k: v + [None] * (max_length - len(v)) for k, v in result_dict.items()}
    return pd.DataFrame(padded_dict)


def export_sequence_list(merged_df, group_data):
    """Export sequence list as CSV bytes."""
    result_df = extract_sequences_by_name(merged_df, group_data)
    if result_df.empty:
        return None
    return result_df.to_csv(index=False).encode('utf-8')


# ---------------------------------------------------------------------------
# Export: Summed peptide results (Excel)
# ---------------------------------------------------------------------------

def summed_peptide_results(merged_df, group_data):
    """Calculate summed peptide results per group."""
    group_data_clean, df = _convert_group_data_dict(merged_df, group_data)
    total_peptide_results_dict = {}
    filtered_df = df.copy()

    for _, value in group_data_clean.items():
        grouping_variable = value['grouping_variable']
        abundance_columns = value['abundance_columns']
        valid_abundance_cols = abundance_columns

        if not valid_abundance_cols:
            continue

        temp_df = filtered_df[['Unique Peptide ID'] + valid_abundance_cols].copy()
        for col in valid_abundance_cols:
            temp_df[col] = pd.to_numeric(temp_df[col], errors='coerce')

        valid_data_mask = (
            temp_df[valid_abundance_cols].notna().any(axis=1) &
            (temp_df[valid_abundance_cols] != 0).any(axis=1) &
            temp_df['Unique Peptide ID'].notna()
        )
        temp_df = temp_df[valid_data_mask]

        if temp_df.empty:
            total_peptide_results_dict[grouping_variable] = {
                'unique_peptides': 0,
                'total_Abundance': 0,
                'total_sem': 0,
                'abundance_sem': 0,
                'count_sem': 0,
                'replicate_data': {
                    'abundance_columns': valid_abundance_cols,
                    'replicate_counts': [0] * len(valid_abundance_cols),
                    'replicate_abundances': [0] * len(valid_abundance_cols)
                }
            }
            continue

        replicate_counts = []
        for col in valid_abundance_cols:
            count = temp_df[temp_df[col].notna() & (temp_df[col] != 0)]['Unique Peptide ID'].nunique()
            replicate_counts.append(count)

        if len(replicate_counts) > 1:
            count_sem = np.std(replicate_counts, ddof=1) / np.sqrt(len(replicate_counts))
        else:
            count_sem = 0

        abundances = temp_df[valid_abundance_cols].values.astype(float)
        peptide_means = np.nanmean(abundances, axis=1)
        total_abundance = np.nansum(peptide_means)

        peptide_sems = np.nanstd(abundances, axis=1) / np.sqrt(abundances.shape[1])
        total_sem = np.sqrt(np.nansum(peptide_sems ** 2))

        all_unique_peptides = temp_df[
            (temp_df[valid_abundance_cols] > 0).any(axis=1)
        ]['Unique Peptide ID'].nunique()

        total_peptide_results_dict[grouping_variable] = {
            'unique_peptides': all_unique_peptides,
            'total_Abundance': total_abundance,
            'total_sem': total_sem,
            'abundance_sem': total_sem,
            'count_sem': count_sem,
            'replicate_data': {
                'abundance_columns': valid_abundance_cols,
                'replicate_counts': replicate_counts,
                'replicate_abundances': [
                    temp_df[col].replace(0, np.nan).sum() for col in valid_abundance_cols
                ]
            }
        }

    return total_peptide_results_dict


def export_summed_peptide_data(merged_df, group_data):
    """Export summed peptide results as Excel bytes."""
    data = summed_peptide_results(merged_df, group_data)
    if not data:
        return None

    summary_data = []
    for group, values in data.items():
        summary_data.append({
            'Group': group,
            'Total_Abundance': values['total_Abundance'],
            'Abundance_SEM': values['abundance_sem'],
            'Unique_Peptides': values['unique_peptides'],
            'Count_SEM': values['count_sem']
        })
    summary_df = pd.DataFrame(summary_data)

    replicate_data = []
    for group, values in data.items():
        replicate_info = values['replicate_data']
        abundance_columns = replicate_info['abundance_columns']
        replicate_abundances = replicate_info['replicate_abundances']
        replicate_counts = replicate_info['replicate_counts']

        for i, replicate_name in enumerate(abundance_columns):
            if i < len(replicate_abundances) and i < len(replicate_counts):
                replicate_data.append({
                    'Group': group,
                    'Replicate': replicate_name,
                    'Total_Abundance': replicate_abundances[i],
                    'Unique_Peptides': replicate_counts[i]
                })
            else:
                replicate_data.append({
                    'Group': group,
                    'Replicate': replicate_name,
                    'Total_Abundance': 0,
                    'Unique_Peptides': 0
                })

    replicate_df = pd.DataFrame(replicate_data)

    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
        replicate_df.to_excel(writer, sheet_name='Replicate Details', index=False)

        for sheet in writer.sheets.values():
            for column_cells in sheet.columns:
                column_cells = [cell for cell in column_cells if cell.value is not None]
                if column_cells:
                    max_length = max(len(str(cell.value)) for cell in column_cells)
                    sheet.column_dimensions[column_cells[0].column_letter].width = max_length + 2
            _apply_scientific_notation(sheet)

    return output.getvalue()


# ---------------------------------------------------------------------------
# Export: Protein analysis (CSV)
# ---------------------------------------------------------------------------

def _fetch_protein_names(accession_str, protein_dict):
    """Fetch protein names from the protein dictionary."""
    names = []
    if '; ' in accession_str:
        accessions = accession_str.split('; ')
    elif '/ ' in accession_str:
        accessions = accession_str.split('/ ')
    else:
        accessions = [accession_str]

    for acc in accessions:
        if acc in protein_dict:
            names.append(protein_dict[acc]['name'])
        else:
            names.append(acc)
    return names


def export_protein_data(merged_df, group_data, protein_dict):
    """
    Export protein analysis results as CSV bytes.

    Returns:
        tuple: (csv_bytes, warnings_list)
    """
    if merged_df is None:
        return None, []

    df = merged_df.copy()
    abundance_cols, selected_groups = _get_abundance_cols(group_data)

    df['Total_Abundance'] = df[abundance_cols].sum(axis=1).astype(int)
    result_df = df[['Unique Peptide ID', 'Total_Abundance']]
    result_df = result_df[result_df['Total_Abundance'] == 0]
    all_zero_list = list(result_df['Unique Peptide ID'])
    peptides_df = df[~df['Unique Peptide ID'].isin(all_zero_list)].copy()

    def normalize_protein(x):
        if isinstance(x, list):
            return ';'.join([str(p) for p in x if p and str(p) != 'nan'])
        elif pd.isna(x):
            return 'Unknown'
        return str(x)

    peptides_df['Protein'] = peptides_df['Protein'].apply(normalize_protein)

    proteins_df = peptides_df.groupby('Protein', as_index=False)[abundance_cols].sum()

    for col in abundance_cols:
        col_sum = proteins_df[col].sum()
        if col_sum > 0:
            proteins_df[f'Rel_{col}'] = (proteins_df[col] / col_sum) * 100
        else:
            proteins_df[f'Rel_{col}'] = 0

    # Add protein descriptions
    name_list = []
    for _, row in proteins_df.iterrows():
        if isinstance(row['Protein'], list):
            strrow = [str(p) for p in row['Protein']]
            named_combo = _fetch_protein_names('; '.join(strrow), protein_dict)
        elif ',' in str(row['Protein']):
            strrow = str(row['Protein']).split(',')
            named_combo = _fetch_protein_names('; '.join(strrow), protein_dict)
        else:
            named_combo = _fetch_protein_names(str(row['Protein']), protein_dict)
        name_list.append(named_combo)

    proteins_df['Description'] = name_list
    proteins_df['Description'] = proteins_df['Description'].astype(str).str.replace(
        r"['\['\]]", "", regex=True
    )

    total_sum = proteins_df[abundance_cols].sum().sum()
    row_sums = proteins_df[abundance_cols].sum(axis=1)
    proteins_df['avg_abundance_all'] = (row_sums / total_sum * 100).round(2) if total_sum > 0 else 0
    proteins_df = proteins_df.sort_values('avg_abundance_all', ascending=False)

    # Track peptide-to-protein mappings
    protein_to_peptides = defaultdict(set)
    protein_to_group_peptides = defaultdict(lambda: defaultdict(set))
    accession_to_idx = {}
    accession_to_description = {}
    for idx, row in proteins_df.iterrows():
        if 'Protein' in row and pd.notna(row['Protein']):
            prot = str(row['Protein'])
            accession_to_idx[prot] = idx
            accession_to_description[prot] = row['Description']
            # For combined proteins (e.g., 'P12345;P67890'), also index each
            # individual accession so the counting loop can find them.
            if ';' in prot:
                for single_acc in [a.strip() for a in prot.split(';') if a.strip()]:
                    if single_acc not in accession_to_idx:
                        accession_to_idx[single_acc] = idx
                        accession_to_description[single_acc] = row['Description']

    warnings_list = []

    for group in selected_groups:
        count_col = f'Count_{group}'
        rel_count_col = f'Rel_Count_{group}'
        proteins_df[count_col] = 0.0
        proteins_df[rel_count_col] = 0.0

        group_peptides = peptides_df[peptides_df[f'Avg_{group}'] > 0]
        peptide_count_total = len(group_peptides)
        counted_peptides = set()

        for _, peptide in group_peptides.iterrows():
            if 'Protein' not in peptide or pd.isna(peptide['Protein']):
                continue
            peptide_id = peptide.get('Unique Peptide ID', None)
            if peptide_id is None or pd.isna(peptide_id):
                continue
            if peptide_id in counted_peptides:
                continue

            accession = peptide['Protein']
            found_match = False

            if ';' in str(accession) or '/' in str(accession):
                sep = ';' if ';' in str(accession) else '/'
                accessions = [acc.strip() for acc in str(accession).split(sep) if acc.strip()]
                for acc in accessions:
                    if acc in accession_to_idx:
                        idx = accession_to_idx[acc]
                        proteins_df.at[idx, count_col] += 1
                        protein_desc = accession_to_description.get(acc, acc)
                        protein_to_peptides[protein_desc].add(peptide_id)
                        protein_to_group_peptides[protein_desc][group].add(peptide_id)
                        counted_peptides.add(peptide_id)
                        found_match = True
                        break
            else:
                if accession in accession_to_idx:
                    idx = accession_to_idx[accession]
                    proteins_df.at[idx, count_col] += 1
                    protein_desc = accession_to_description.get(accession, accession)
                    protein_to_peptides[protein_desc].add(peptide_id)
                    protein_to_group_peptides[protein_desc][group].add(peptide_id)
                    counted_peptides.add(peptide_id)
                    found_match = True

        if peptide_count_total > 0:
            for idx in proteins_df.index:
                protein_count = float(proteins_df.at[idx, count_col])
                proteins_df.at[idx, rel_count_col] = (protein_count / peptide_count_total) * 100

        warnings_list.append({
            'group': group,
            'total': len(group_peptides),
            'counted': len(counted_peptides)
        })

    # Save Protein ID before dropping the Protein column so it appears in all sheets.
    if 'Protein' in proteins_df.columns:
        proteins_df.insert(0, 'Protein ID', proteins_df['Protein'])
        proteins_df = proteins_df.drop(columns=['Protein'])

    # Split into four sheets by column type
    all_cols = list(proteins_df.columns)
    id_cols   = ['Protein ID']   if 'Protein ID'   in all_cols else []
    desc_cols = ['Description']  if 'Description'  in all_cols else []
    meta_cols = id_cols + desc_cols

    abs_abs_cols  = meta_cols + [c for c in all_cols if c.startswith('Avg_')] + \
                    (['avg_abundance_all'] if 'avg_abundance_all' in all_cols else [])
    abs_rel_cols  = meta_cols + [c for c in all_cols if c.startswith('Rel_Avg_')]
    cnt_abs_cols  = meta_cols + [c for c in all_cols if c.startswith('Count_')]
    cnt_rel_cols  = meta_cols + [c for c in all_cols if c.startswith('Rel_Count_')]

    # Round relative sheets to 2 dp; ensure absolute count columns are integers.
    rel_num_cols_abs = [c for c in abs_rel_cols if c not in meta_cols]
    rel_num_cols_cnt = [c for c in cnt_rel_cols if c not in meta_cols]
    cnt_num_cols     = [c for c in cnt_abs_cols  if c not in meta_cols]

    abs_rel_df = proteins_df[abs_rel_cols].copy()
    abs_rel_df[rel_num_cols_abs] = abs_rel_df[rel_num_cols_abs].round(2)

    cnt_rel_df = proteins_df[cnt_rel_cols].copy()
    cnt_rel_df[rel_num_cols_cnt] = cnt_rel_df[rel_num_cols_cnt].round(2)

    cnt_abs_df = proteins_df[cnt_abs_cols].copy()
    cnt_abs_df[cnt_num_cols] = cnt_abs_df[cnt_num_cols].round(0).astype('Int64')

    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine='openpyxl') as writer:
        proteins_df[abs_abs_cols].to_excel(writer, sheet_name='Abundance Absolute', index=False)
        abs_rel_df.to_excel(writer, sheet_name='Abundance Relative', index=False)
        cnt_abs_df.to_excel(writer, sheet_name='Count Absolute',      index=False)
        cnt_rel_df.to_excel(writer, sheet_name='Count Relative',       index=False)
        for ws in writer.sheets.values():
            _apply_scientific_notation(ws)

    return buf.getvalue(), warnings_list


# ---------------------------------------------------------------------------
# Export: Summed function data (Excel)
# ---------------------------------------------------------------------------

def _process_bioactive_data(merged_df, group_data):
    """Process bioactive peptide data for functional analysis."""
    if merged_df is None:
        return None

    df = merged_df.copy()
    abundance_cols = []
    selected_groups = []
    unique_function_abundance = {}

    for _, value in group_data.items():
        grouping_variable = value['grouping_variable']
        abundance_columns = value['abundance_columns']
        selected_groups.append(grouping_variable)
        if grouping_variable != abundance_columns:
            abundance_cols.append(f'Avg_{grouping_variable}')
        else:
            abundance_cols.append(abundance_columns)

        if 'function' not in df.columns:
            return None

        for column in abundance_cols:
            temp_df = df[['Unique Peptide ID', 'function', column]].copy()
            temp_df = temp_df[
                (temp_df[column] != 0) &
                temp_df[column].notna() &
                temp_df['function'].notna()
            ]
            if temp_df.empty:
                continue

            temp_df.loc[:, 'function'] = temp_df['function'].fillna('').str.split(';')
            exploded_df = temp_df.explode('function')
            exploded_df.loc[:, 'function'] = exploded_df['function'].str.strip()
            exploded_df = exploded_df[exploded_df['function'] != '']

            if not exploded_df.empty:
                function_grouped = exploded_df.groupby('function')[column].sum()
                unique_function_abundance[grouping_variable] = function_grouped.to_dict()

    return unique_function_abundance, abundance_cols


def _process_functional_peptide_export_data(merged_df, group_data):
    """Process data for functional peptide export into Excel format."""
    result = _process_bioactive_data(merged_df, group_data)
    if result is None:
        return None
    unique_function_abundance, abundance_cols = result
    if not unique_function_abundance:
        return None

    groups = list(unique_function_abundance.keys())
    df = merged_df.copy()

    summed_function_count = {}
    unique_function_counts = {}
    summed_function_abundance = {}

    for group, abundance_column in zip(groups, abundance_cols):
        if abundance_column not in df.columns:
            continue

        temp_df = df[['Unique Peptide ID', 'function', abundance_column]].copy()
        temp_df = temp_df[
            (temp_df[abundance_column] != 0) &
            temp_df[abundance_column].notna() &
            temp_df['function'].notna()
        ]

        filtered_df = temp_df.drop_duplicates(subset='Unique Peptide ID')
        unique_peptide_count = filtered_df['Unique Peptide ID'].nunique()
        total_sum = filtered_df[abundance_column].sum()

        summed_function_abundance[group] = total_sum
        summed_function_count[group] = unique_peptide_count

        filtered_df = filtered_df.copy()
        filtered_df.loc[:, 'function'] = filtered_df['function'].fillna('').str.split(';')
        exploded_df = filtered_df.explode('function')
        exploded_df.loc[:, 'function'] = exploded_df['function'].str.strip()
        exploded_df = exploded_df[exploded_df['function'] != '']

        if not exploded_df.empty:
            function_counts = exploded_df['function'].value_counts().to_dict()
            unique_function_counts[group] = function_counts

    peptide_count_df = pd.DataFrame.from_dict(
        summed_function_count, orient='index', columns=['Counts of peptides']
    )
    function_count_df = pd.DataFrame.from_dict(
        unique_function_counts, orient='index'
    ).fillna(0).astype(int)
    combined_count_df = pd.concat([peptide_count_df, function_count_df], axis=1).T

    peptide_abundance_df = pd.DataFrame.from_dict(
        summed_function_abundance, orient='index', columns=['Summed Abundance']
    )
    function_abundance_df = pd.DataFrame.from_dict(
        unique_function_abundance, orient='index'
    ).fillna(0)
    combined_abundance_df = pd.concat(
        [peptide_abundance_df, function_abundance_df], axis=1
    ).T

    combined_df = pd.DataFrame(
        index=combined_abundance_df.index,
        columns=combined_abundance_df.columns
    )
    for col in combined_abundance_df.columns:
        for idx in combined_abundance_df.index:
            abundance = combined_abundance_df.loc[idx, col]
            count = (combined_count_df.loc['Counts of peptides', col]
                     if idx == 'Summed Abundance'
                     else combined_count_df.loc[idx, col])
            combined_df.loc[idx, col] = (
                "-" if (abundance == 0 and count == 0)
                else f"{abundance:.2e} ({round(count)})"
            )
    combined_df.rename(index={'Summed Abundance': 'Total'}, inplace=True)

    return combined_df, combined_count_df, combined_abundance_df


def export_summed_function_data(merged_df, group_data):
    """Export summed function data as Excel bytes with absolute and relative sheets."""
    result = _process_functional_peptide_export_data(merged_df, group_data)
    if result is None:
        return None

    combined_df, combined_count_df, combined_abundance_df = result

    # --- Abundance Relative ---
    # Each function's % = function_abundance / sum_of_all_function_abundances * 100.
    # Using the function-row sum (not the unique-peptide total) ensures percentages
    # sum to exactly 100%, even when peptides carry multiple functions.
    abundance_rel_df = combined_abundance_df.copy().astype(float)
    _abs_total_label = 'Summed Abundance'
    _abs_func_rows = abundance_rel_df.index[abundance_rel_df.index != _abs_total_label]
    if len(_abs_func_rows) > 0:
        func_totals = abundance_rel_df.loc[_abs_func_rows].sum(axis=0)
        for col in abundance_rel_df.columns:
            denom = func_totals[col] if func_totals[col] != 0 else np.nan
            abundance_rel_df.loc[_abs_func_rows, col] = (
                abundance_rel_df.loc[_abs_func_rows, col] / denom * 100
            )
    if _abs_total_label in abundance_rel_df.index:
        abundance_rel_df.loc[_abs_total_label] = 100.0

    # --- Count Relative ---
    # Each function's % = function_count / sum_of_all_function_counts * 100.
    count_rel_df = combined_count_df.copy().astype(float)
    _cnt_total_label = 'Counts of peptides'
    _cnt_func_rows = count_rel_df.index[count_rel_df.index != _cnt_total_label]
    if len(_cnt_func_rows) > 0:
        func_totals = count_rel_df.loc[_cnt_func_rows].sum(axis=0)
        for col in count_rel_df.columns:
            denom = func_totals[col] if func_totals[col] != 0 else np.nan
            count_rel_df.loc[_cnt_func_rows, col] = (
                count_rel_df.loc[_cnt_func_rows, col] / denom * 100
            )
    if _cnt_total_label in count_rel_df.index:
        count_rel_df.loc[_cnt_total_label] = 100.0

    # Round relative sheets to 2 dp; count absolute should be whole numbers.
    abundance_rel_rounded = abundance_rel_df.round(2)
    count_rel_rounded = count_rel_df.round(2)
    count_abs_rounded = combined_count_df.copy()
    count_abs_rounded = count_abs_rounded.round(0)

    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        combined_df.to_excel(writer, sheet_name='Combined', index=True)
        combined_abundance_df.to_excel(writer, sheet_name='Abundance Absolute', index=True)
        abundance_rel_rounded.to_excel(writer, sheet_name='Abundance Relative', index=True)
        count_abs_rounded.to_excel(writer, sheet_name='Count Absolute', index=True)
        count_rel_rounded.to_excel(writer, sheet_name='Count Relative', index=True)
        for ws in writer.sheets.values():
            _apply_scientific_notation(ws)

    return output.getvalue()


# ---------------------------------------------------------------------------
# Export: Technical replicate (within-biological-replicate) correlation (Excel)
# ---------------------------------------------------------------------------

def export_tech_rep_correlation(pd_results_df, tech_dup_mapping,
                                correlation_type='Pearson', log_transform=True,
                                log_base=10):
    """
    Calculate and export pairwise correlations between technical replicates.

    For each biological replicate produced by collapsing technical duplicates,
    computes pairwise Pearson/Spearman correlations between the original
    technical replicate columns in the pre-collapsed peptidomic data.

    Args:
        pd_results_df: original peptidomic DataFrame (pre-collapse) containing
                       the raw technical replicate columns
        tech_dup_mapping: {biological_replicate_name: [tech_rep_col1, tech_rep_col2, ...]}
        correlation_type: 'Pearson' or 'Spearman'
        log_transform: apply log10 transform before correlating

    Returns:
        Excel bytes or None if no valid pairs exist
    """
    if not tech_dup_mapping or pd_results_df is None:
        return None

    df = pd_results_df.copy()
    within_rep_correlations = {}

    for bio_rep, tech_rep_cols in tech_dup_mapping.items():
        valid_cols = [c for c in tech_rep_cols if c in df.columns]
        if len(valid_cols) < 2:
            continue

        data = df[valid_cols].copy()
        for col in valid_cols:
            data[col] = pd.to_numeric(data[col], errors='coerce')
        data = data[data.gt(0).all(axis=1)]

        if len(data) < 2:
            continue

        data = prepare_data(data, log_transform, log_base)
        pairs, correlations, p_values = [], [], []

        for i in range(len(valid_cols)):
            for j in range(i):
                col1, col2 = valid_cols[j], valid_cols[i]
                mask = data[[col1, col2]].notnull().all(axis=1) & \
                       (data[[col1, col2]] > 0).all(axis=1)
                x, y = data.loc[mask, col1], data.loc[mask, col2]
                if len(x) > 1:
                    corr, pval = calculate_correlation(x, y, correlation_type)
                    correlations.append(round(corr, 3))
                    p_values.append(pval)
                else:
                    correlations.append(np.nan)
                    p_values.append(np.nan)
                pairs.append(f"{col1} vs {col2}")

        within_rep_correlations[bio_rep] = pd.DataFrame({
            'Pair': pairs,
            'Correlation': correlations,
            'P-Value': p_values,
        })

    if not within_rep_correlations:
        return None

    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        # All Correlations sheet
        combined_rows = []
        for bio_rep, df_corr in within_rep_correlations.items():
            temp = df_corr.copy()
            temp.insert(0, 'Biological Replicate', bio_rep)
            combined_rows.append(temp)
        pd.concat(combined_rows, ignore_index=True).to_excel(
            writer, sheet_name='All Correlations', index=False
        )

        # Summary sheet
        summary_list = []
        for bio_rep, df_corr in within_rep_correlations.items():
            valid_corrs = df_corr['Correlation'].dropna()
            summary_list.append({
                'Biological Replicate': bio_rep,
                'Technical Replicates': ', '.join(tech_dup_mapping[bio_rep]),
                'N Technical Replicates': len(tech_dup_mapping[bio_rep]),
                'Avg Correlation': round(valid_corrs.mean(), 3) if not valid_corrs.empty else np.nan,
                'Min Correlation': round(valid_corrs.min(), 3) if not valid_corrs.empty else np.nan,
                'Max Correlation': round(valid_corrs.max(), 3) if not valid_corrs.empty else np.nan,
                'N Comparisons': len(valid_corrs),
            })
        _write_summary_with_transform(
            writer, pd.DataFrame(summary_list),
            log_label(log_transform, log_base), correlation_type)

        # One sheet per biological replicate
        for bio_rep, df_corr in within_rep_correlations.items():
            sheet_name = bio_rep[:31]  # Excel sheet name limit
            df_corr.to_excel(writer, sheet_name=sheet_name, index=False)

    return buffer.getvalue()


# ---------------------------------------------------------------------------
# Export: Group (sample-to-sample) correlation (Excel)
# ---------------------------------------------------------------------------

def export_group_correlation(merged_df, group_data, correlation_type='Pearson',
                             log_transform=True, log_base=10):
    """
    Calculate and export cross-group correlation analysis to Excel bytes.
    """
    df = merged_df.copy()
    correlation_results = []
    avg_columns = {
        group_info['grouping_variable']: f"Avg_{group_info['grouping_variable']}"
        for group_info in group_data.values()
        if f"Avg_{group_info['grouping_variable']}" in df.columns
    }

    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        for (group1, col1), (group2, col2) in combinations(avg_columns.items(), 2):
            mask = (df[col1] > 0) & (df[col2] > 0)
            if mask.sum() > 1:
                values1 = prepare_data(df.loc[mask, col1], log_transform, log_base)
                values2 = prepare_data(df.loc[mask, col2], log_transform, log_base)
                correlation, p_value = calculate_correlation(values1, values2, correlation_type)
                correlation_results.append({
                    'Group 1': group1,
                    'Group 2': group2,
                    'Correlation': round(correlation, 3),
                    'P-Value': p_value,
                    'Number of Peptides': mask.sum()
                })

        if correlation_results:
            cross_group_df = pd.DataFrame(correlation_results)
            cross_group_df.to_excel(writer, sheet_name='Cross-Group Correlations', index=False)

            summary_stats = {
                'Average Correlation': round(np.mean([r['Correlation'] for r in correlation_results]), 3),
                'Min Correlation': round(min([r['Correlation'] for r in correlation_results]), 3),
                'Max Correlation': round(max([r['Correlation'] for r in correlation_results]), 3),
                'Total Comparisons': len(correlation_results)
            }
            _write_summary_with_transform(
                writer, pd.DataFrame([summary_stats]),
                log_label(log_transform, log_base), correlation_type)
        else:
            pd.DataFrame({'Message': ['No valid group pairs for correlation']}).to_excel(
                writer, sheet_name='Summary', index=False
            )

    return buffer.getvalue()


# ---------------------------------------------------------------------------
# Export: Replicate (within-group) correlation (Excel)
# ---------------------------------------------------------------------------

def export_replicate_correlation(merged_df, group_data, correlation_type='Pearson',
                                  log_transform=True, log_base=10):
    """
    Calculate and export within-group replicate correlation analysis to Excel bytes.
    """
    group_data_clean, df = _convert_group_data_dict(merged_df, group_data)
    within_group_correlations = {}

    for key, value in group_data_clean.items():
        grouping_variable = value['grouping_variable']
        abundance_columns = value['abundance_columns']

        data = df[abundance_columns].copy()
        data = data[data.gt(0).all(axis=1)]

        if len(data) > 1:
            data = prepare_data(data, log_transform, log_base)
            pairs = []
            correlations = []
            p_values = []

            n_cols = len(abundance_columns)
            for i in range(n_cols):
                for j in range(i):
                    col1 = abundance_columns[j]
                    col2 = abundance_columns[i]
                    pair_name = f"{col1} vs {col2}"
                    mask = data[[col1, col2]].notnull().all(axis=1) & (data[[col1, col2]] > 0).all(axis=1)
                    x = data.loc[mask, col1]
                    y = data.loc[mask, col2]
                    if len(x) > 1 and len(y) > 1:
                        corr, pval = calculate_correlation(x, y, correlation_type)
                        correlations.append(round(corr, 3))
                        p_values.append(pval)
                    else:
                        correlations.append(np.nan)
                        p_values.append(np.nan)
                    pairs.append(pair_name)

            within_group_correlations[grouping_variable] = pd.DataFrame({
                'Pair': pairs,
                'Correlation': correlations,
                'P-Value': p_values
            })

    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
        if within_group_correlations:
            combined_corrs = []
            for group, df_corr in within_group_correlations.items():
                temp = df_corr.copy()
                temp['Group'] = group
                combined_corrs.append(temp)
            combined_correlation_df = pd.concat(combined_corrs, ignore_index=True)

            summary_list = []
            for group, df_corr in within_group_correlations.items():
                valid_corrs = df_corr['Correlation'].dropna()
                if not valid_corrs.empty:
                    summary_list.append({
                        'Group': group,
                        'Min Correlation': round(valid_corrs.min(), 3),
                        'Max Correlation': round(valid_corrs.max(), 3),
                        'Average Correlation': round(valid_corrs.mean(), 3),
                        'Total Comparisons': len(valid_corrs)
                    })
                else:
                    summary_list.append({
                        'Group': group,
                        'Min Correlation': np.nan,
                        'Max Correlation': np.nan,
                        'Average Correlation': np.nan,
                        'Total Comparisons': 0
                    })
            summary_df = pd.DataFrame(summary_list)

            combined_correlation_df.to_excel(writer, sheet_name='All Correlations', index=False)
            _write_summary_with_transform(
                writer, summary_df, log_label(log_transform, log_base),
                correlation_type)

            for key, value in group_data_clean.items():
                grouping_variable = value['grouping_variable']
                if grouping_variable in within_group_correlations:
                    group_df = within_group_correlations[grouping_variable]
                    sheet_name = grouping_variable[:31]  # Excel sheet name limit
                    group_df.to_excel(writer, sheet_name=sheet_name, index=False)
        else:
            pd.DataFrame({'Message': ['No valid replicate pairs for correlation']}).to_excel(
                writer, sheet_name='Summary', index=False
            )

    return buffer.getvalue()
