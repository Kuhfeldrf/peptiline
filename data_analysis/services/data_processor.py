"""
Data processing service for the Data Analysis web app.
Extracted from data_analysis.ipynb DataTransformation and DataHandler classes.
"""
import re
import io
import traceback

import numpy as np
import pandas as pd


# Protein-name display cleaning, shared verbatim with the Heatmap app
# (heatmap_viz/services/data_processor._clean_protein_name) so the "Strip protein
# name" toggle behaves identically across both applications.
#
# Two kinds of clutter are removed:
#   * Trailing UniProt FASTA metadata — "Beta-casein OS=Bos taurus GN=CSN2 …".
#   * A leading UniProt entry-name token — "LACB_BOVIN Beta-lactoglobulin",
#     "B4GT1_BOVIN Beta-1,4-galactosyltransferase 1" (the "CAS_bovine" leader).
# Entry names are always upper-case ID_SPECIES, so the leader regex is upper-case
# only and requires a descriptive name after it — a bare "LACB_BOVIN" (no space)
# and lower-cased names are left untouched, so ordinary names never get clipped.
_FASTA_META_RE = re.compile(r'\s+(?:OS|OX|GN|PE|SV)=')
_ENTRY_NAME_LEADER_RE = re.compile(r'^[A-Z0-9]+_[A-Z0-9]+\s+(?=\S)')


def _clean_protein_name(name) -> str:
    """Strip trailing UniProt FASTA metadata and a leading entry-name token."""
    s = _FASTA_META_RE.split(str(name), 1)[0].strip()
    return _ENTRY_NAME_LEADER_RE.sub('', s, count=1).strip()


_PROTEIN_DELIM_RE = re.compile(r'\s*;\s*|\s*/\s*|\s*,\s*')


def _split_protein_ids(value) -> list:
    """Split one 'Protein' cell into individual protein IDs.

    Cells are either a delimited string (';', '/', ',') from a directly
    uploaded CSV, or a genuine Python list object (from a pickled/transferred
    dataframe — Data Transformation's extract_protein_id() stores a real list
    for multi-protein peptides; see extract_protein_dict() above for why a
    bare str() would mangle that case).
    """
    if isinstance(value, list):
        return [str(p).strip() for p in value if str(p).strip()]
    if pd.isna(value):
        return []
    return [p for p in _PROTEIN_DELIM_RE.split(str(value)) if p]


def _compute_protein_abundance(df: pd.DataFrame, avg_columns: list) -> dict:
    """
    Total abundance per protein ID, vectorized.

    Mirrors the old row-by-row version (each peptide's abundance split evenly
    across the protein IDs on its row) without an `iterrows()` pass: the
    per-row abundance sum uses pandas' vectorized `.sum(axis=1)`, and the
    per-protein accumulation uses `explode()` + `groupby().sum()` instead of a
    Python dict built up one row at a time. On wide/tall uploads this is the
    difference between a sub-second call and one that can run for minutes.
    """
    if 'Protein' not in df.columns or not avg_columns:
        return {}

    total_ab = df[avg_columns].apply(pd.to_numeric, errors='coerce').fillna(0).sum(axis=1)
    parts = df['Protein'].map(_split_protein_ids)
    n_parts = parts.map(len)

    has_parts = n_parts > 0
    if not has_parts.any():
        return {}

    per_protein = pd.Series(0.0, index=df.index)
    per_protein[has_parts] = total_ab[has_parts] / n_parts[has_parts]

    exploded = pd.DataFrame({'pid': parts, 'val': per_protein}).explode('pid')
    exploded = exploded.dropna(subset=['pid'])
    if exploded.empty:
        return {}
    return exploded.groupby('pid')['val'].sum().to_dict()


_BROADER_FUNCTIONS = {'Functional Peptides', 'Non-Functional Peptides', 'All Functional Peptides', 'Other Functions'}


def _compute_function_counts(df: pd.DataFrame) -> dict:
    """
    Peptide count per individual function label.

    NOTE: measured, this is deliberately *not* an explode()/value_counts()
    pipeline. That was tried first and benchmarked slower than this plain
    loop at every size from 50k to 500k rows (pandas' per-call overhead on
    split/explode/dropna/strip/isin outweighs what a tight Python loop over
    already-split strings costs). Unlike the protein-abundance and
    protein-dict cases above, this one was never an iterrows() bottleneck --
    dropna() already does the heavy lifting, and iterating the leftover
    Series of strings is cheap. Kept as a shared function purely to
    de-duplicate the two call sites (get_selector_options and
    DataAnalysisState._resolve_all_functions).
    """
    if 'function' not in df.columns:
        return {}
    function_totals: dict = {}
    for func_str in df['function'].dropna():
        if isinstance(func_str, str):
            for func in [f.strip() for f in func_str.split(';')]:
                if func and func not in _BROADER_FUNCTIONS:
                    function_totals[func] = function_totals.get(func, 0) + 1
    return function_totals


# ---------------------------------------------------------------------------
# Selection sentinels + plot-filter derivation
# ---------------------------------------------------------------------------
# The UI used to carry an explicit "Plot Filter" dropdown that decided which of
# the protein / function selectors actually applied; the selectors themselves
# were greyed out until the dropdown said otherwise, which users consistently
# missed. The dropdown is gone: what the user selects IS the filter. Each
# selector carries an explicit "None" sentinel at the top so that "no filter on
# this dimension" stays distinguishable from "every protein/function", which are
# genuinely different plots (the latter breaks the data out per item).

NO_PROTEIN_FILTER = 'No Protein Filter'
ALL_PROTEINS = 'All Proteins (No Filter)'
NO_FUNCTION_FILTER = 'No Function Filter'
ALL_FUNCTIONAL = 'All Functional Peptides'
NON_FUNCTIONAL = 'Non-Functional Peptides'


def _selection_is_active(raw, none_sentinel) -> bool:
    """True when a selector actually restricts/breaks out the plot."""
    return bool([v for v in (raw or []) if v and v != none_sentinel])


def derive_plot_filter(selected_proteins_raw, selected_functions_raw) -> str:
    """Infer the legacy plot_filter mode from the two selections.

    Kept as a string because the plotting layer branches on these five modes
    throughout; only its *source* changed (selections instead of a dropdown).
    Selecting both "All Functional Peptides" and "Non-Functional Peptides" is
    how the old "Functional vs Non-Functional" mode is now requested.
    """
    prot = _selection_is_active(selected_proteins_raw, NO_PROTEIN_FILTER)
    func = _selection_is_active(selected_functions_raw, NO_FUNCTION_FILTER)

    fn_set = set(selected_functions_raw or [])
    if func and ALL_FUNCTIONAL in fn_set and NON_FUNCTIONAL in fn_set:
        return 'Functional vs Non-Functional Peptides'

    if prot and func:
        return 'Both'
    if prot:
        return 'Selected Protein(s)'
    if func:
        return 'Selected Function(s)'
    return 'No Filter'


# ---------------------------------------------------------------------------
# File loading
# ---------------------------------------------------------------------------

def load_file(file_obj, filename: str):
    """Load a CSV/TSV/XLSX/Parquet file, return (df, error_message)."""
    name_lower = filename.lower()
    try:
        if name_lower.endswith('.parquet'):
            df = pd.read_parquet(file_obj)
        elif name_lower.endswith('.xlsx'):
            df = pd.read_excel(file_obj)
        elif name_lower.endswith('.tsv') or name_lower.endswith('.txt'):
            df = pd.read_csv(file_obj, sep='\t', low_memory=False)
        else:
            df = pd.read_csv(file_obj, low_memory=False)
        df.columns = df.columns.str.strip()
        # Excel-sourced exports can carry thousands of fully-blank trailing rows
        # (all-comma/all-empty lines past the real data). Their columns parse as
        # NaN/float while the real rows parse as str, which is what triggers
        # pandas' "mixed types" DtypeWarning — and every downstream iterrows()/
        # apply() pass over the dataframe then wastes time on rows with no data.
        # Same fix already used by the Data Transformation loader (data_loader.py).
        df = df.dropna(how='all')
        df = df[~df.astype(str).apply(lambda row: row.str.strip().eq('').all(), axis=1)]
        return df, None
    except Exception as exc:
        return None, str(exc)


# ---------------------------------------------------------------------------
# Group data extraction
# ---------------------------------------------------------------------------

def process_group_data(df: pd.DataFrame):
    """
    Extract group data from DataFrame columns.
    Handles both 'Grouped:' pattern and plain 'Avg_*' columns.
    Returns (group_data_dict, renamed_df, warning).
      group_data_dict = {group_name: [col1, col2, ...]}  (list of replicate cols)
                     OR {group_name: 'Avg_GroupName'}    (no-replicate form)
    """
    group_data_dict = {}
    renamed_columns = {}
    warning = None

    grouped_columns = [col for col in df.columns if " 'Grouped:" in str(col)]

    if not grouped_columns:
        avg_columns = [col for col in df.columns if col.startswith('Avg_')]
        if not avg_columns:
            return {}, df, "No group columns (Avg_* or Grouped:) found in file."
        for col in avg_columns:
            group_name = col.replace('Avg_', '')
            group_data_dict[group_name] = col
        warning = "No replicate data found. Some features (error bars, correlation) will be unavailable."
        return group_data_dict, df, warning

    for col in grouped_columns:
        base_col_name = col.split(" 'Grouped:")[0].strip()
        match = re.search(r"\((.*?)\)", col)
        if match:
            groups_str = match.group(1)
            groups = [g.strip() for g in groups_str.split(";")]
            for group in groups:
                if group not in group_data_dict:
                    group_data_dict[group] = []
                group_data_dict[group].append(base_col_name)
            renamed_columns[col] = base_col_name

    df_renamed = df.rename(columns=renamed_columns)
    return group_data_dict, df_renamed, warning


# ---------------------------------------------------------------------------
# Protein info extraction
# ---------------------------------------------------------------------------

def extract_protein_dict(df: pd.DataFrame) -> dict:
    """Build {protein_id: {name, species, description}} from DataFrame.

    Vectorized: called unconditionally on every upload/transfer (same as
    get_selector_options), so this used to be a second full iterrows() pass
    over the whole peptide table. explode() replaces the per-row ID-splitting
    loop, and a single drop_duplicates(keep='first') reproduces the old "first
    occurrence wins" rule for name/species (rows are exploded in original
    order, so the first row a protein ID appears in is still what wins).
    """
    protein_dict: dict = {}
    if 'Protein' not in df.columns:
        return protein_dict

    # Data Transformation's extract_protein_id() stores a genuine Python list
    # (not a joined string) for multi-protein peptides -- _split_protein_ids
    # handles both that and the delimited-string case from a directly
    # uploaded CSV (see its docstring for why str()-ing a list would mangle it).
    parts = df['Protein'].map(_split_protein_ids)
    has_parts = parts.map(len) > 0
    if not has_parts.any():
        return protein_dict

    # NA sentinels reproduce row.get(col, default) when the column is absent
    # entirely: pd.isna() below then falls back exactly as the old code did.
    name_col = df['protein_name'] if 'protein_name' in df.columns else pd.Series(pd.NA, index=df.index)
    species_col = df['protein_species'] if 'protein_species' in df.columns else pd.Series(pd.NA, index=df.index)

    exploded = pd.DataFrame({
        'pid': parts[has_parts],
        'name_value': name_col[has_parts],
        'species_value': species_col[has_parts],
    }).explode('pid')
    exploded = exploded.drop_duplicates(subset='pid', keep='first')

    for pid, name_value, species_value in zip(
        exploded['pid'], exploded['name_value'], exploded['species_value']
    ):
        # Data Transformation's own "cleaning up placeholder values" step
        # (data_combiner.py) replaces the literal string 'Unknown' with
        # pd.NA before the merged dataframe is saved -- so a transferred
        # dataset's unresolved names arrive as pd.NA, not the string
        # 'Unknown'. pd.isna() catches that (and None/NaN) uniformly, which
        # plain string matching cannot. Treat any of those, or the literal
        # string itself (belt-and-suspenders for other sources), as no name
        # at all and fall back to the protein ID -- matches the Heatmap
        # module's convention so an unmapped protein is never shown as a
        # bare, unhelpful placeholder.
        if pd.isna(name_value):
            name = pid
        else:
            name = str(name_value).strip()
            if not name or name.lower() in ('unknown', 'nan', 'none'):
                name = pid
        species = 'Unknown' if pd.isna(species_value) else str(species_value)
        protein_dict[pid] = {'name': name, 'species': species}
    return protein_dict


# ---------------------------------------------------------------------------
# Validate / standardize columns
# ---------------------------------------------------------------------------

def validate_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, list]:
    """Ensure required columns exist. Returns (df, list_of_warnings).

    Matching is case-insensitive (headers are already whitespace-stripped in
    load_file()). This does not extend to the 'Grouped:'/'Avg_*' markers,
    which are generated internally by the app in a fixed case and are never
    user-authored column names.
    """
    warnings = []
    required = ['Unique Peptide ID', 'Protein']
    for col in required:
        if col not in df.columns:
            for existing in df.columns:
                if existing.lower() == col.lower():
                    df = df.rename(columns={existing: col})
                    break
            else:
                if col != 'Protein':
                    warnings.append(f"Required column '{col}' not found.")
    # Ensure protein_name column (may come from different source)
    if 'protein_name' not in df.columns and 'Protein' in df.columns:
        df['protein_name'] = df['Protein']
    return df, warnings


# ---------------------------------------------------------------------------
# Selector options for the UI
# ---------------------------------------------------------------------------

def get_selector_options(merged_df: pd.DataFrame, group_data_dict: dict, protein_dict: dict) -> dict:
    """
    Compute options for the frontend selectors.
    Returns dict with groups, proteins, functions, has_functions.
    """
    groups = list(group_data_dict.keys())
    avg_columns = [col for col in merged_df.columns if col.startswith('Avg_')]

    # Proteins – sorted by total abundance
    protein_abundance = _compute_protein_abundance(merged_df, avg_columns)

    all_proteins_sorted = sorted(
        protein_dict.keys(),
        key=lambda p: protein_abundance.get(p, 0),
        reverse=True,
    )
    protein_options = [
        {'id': pid, 'label': protein_dict[pid].get('name', pid)}
        for pid in all_proteins_sorted
    ]

    # Functions
    has_functions = (
        'function' in merged_df.columns
        and not merged_df['function'].isna().all()
    )
    function_totals = _compute_function_counts(merged_df) if has_functions else {}
    functions = [f for f, _ in sorted(function_totals.items(), key=lambda x: x[1], reverse=True)]

    # Which groups carry replicate-level ('Grouped:') columns. In this module a
    # group's value is a list of replicate column names (Grouped form) or a plain
    # 'Avg_<group>' string (no-replicate form); ≥2 replicates are needed for the
    # SEM / group-comparison statistics. Surfaced so the UI can show the same
    # "replicate data detected / disabled" banner as the Heatmap dashboard.
    var_replicates = {
        g: (isinstance(cols, (list, tuple)) and len(cols) >= 2)
        for g, cols in group_data_dict.items()
    }
    has_replicates = any(var_replicates.values())

    return {
        'groups': groups,
        'proteins': protein_options,
        'protein_ids': all_proteins_sorted,
        'functions': functions,
        'has_functions': has_functions,
        'avg_columns': avg_columns,
        'var_replicates': var_replicates,
        'has_replicates': has_replicates,
    }


# ---------------------------------------------------------------------------
# Colour helpers
# ---------------------------------------------------------------------------

def get_color_sequence(n_colors: int, scheme: str = 'HSV') -> list:
    """Generate n_colors colours from the given plotly/custom scheme."""
    import plotly.express as px

    if n_colors <= 0:
        return []

    single_colors = {
        'red', 'green', 'blue', 'yellow', 'purple', 'orange', 'cyan',
        'magenta', 'pink', 'brown', 'black', 'white', 'gray', 'darkblue',
        'darkgreen', 'darkred', 'darkorange', 'darkpurple', 'lightblue',
        'lightgreen', 'lightred', 'gold', 'silver', 'teal', 'navy', 'maroon',
        'olive', 'lime', 'aqua', 'indigo', 'violet', 'turquoise', 'coral',
        'crimson', 'salmon', 'sienna', 'tan', 'khaki', 'plum', 'orchid',
        'steelblue', 'seagreen', 'mediumpurple', 'slategray',
    }

    try:
        if scheme.startswith('---'):
            scheme = 'HSV'
        elif scheme.lower() in single_colors:
            return [scheme.lower()] * n_colors

        if scheme.lower() in ('rainbow', 'hsv'):
            return [f'hsl({h},70%,60%)' for h in np.linspace(0, 330, n_colors)]

        for palette_group in (
            px.colors.qualitative,
            px.colors.sequential,
            px.colors.diverging,
            px.colors.cyclical,
        ):
            color_sequence = getattr(palette_group, scheme, None)
            if color_sequence:
                if n_colors >= len(color_sequence):
                    indices = np.linspace(0, len(color_sequence) - 1, n_colors)
                    return [color_sequence[int(i)] for i in indices]
                indices = np.linspace(0, len(color_sequence) - 1, n_colors, dtype=int)
                return [color_sequence[i] for i in indices]

        return [f'hsl({h},70%,60%)' for h in np.linspace(0, 330, n_colors)]
    except Exception:
        return [f'hsl({h},70%,60%)' for h in np.linspace(0, 330, n_colors)]


def get_single_color(scheme: str = 'HSV') -> str:
    """Return a single representative colour for the scheme.
    For named single-colour schemes, returns that colour.
    For palette/gradient schemes, returns Carolina Blue (#4B9CD3).
    """
    if scheme.startswith('---'):
        scheme = 'HSV'
    single_colors = {
        'red', 'green', 'blue', 'yellow', 'purple', 'orange', 'cyan',
        'magenta', 'pink', 'brown', 'black', 'white', 'gray', 'darkblue',
        'darkgreen', 'darkred', 'darkorange', 'darkpurple', 'lightblue',
        'lightgreen', 'lightred', 'gold', 'silver', 'teal', 'navy', 'maroon',
        'olive', 'lime', 'aqua', 'indigo', 'violet', 'turquoise', 'coral',
        'crimson', 'salmon', 'sienna', 'tan', 'khaki', 'plum', 'orchid',
        'steelblue', 'seagreen', 'mediumpurple', 'slategray',
    }
    if scheme.lower() in single_colors:
        return scheme.lower()
    return '#4B9CD3'  # Carolina Blue default for palette/gradient/HSV schemes


# ---------------------------------------------------------------------------
# Core: contains_function helper
# ---------------------------------------------------------------------------

def contains_function(func_string, target_function: str) -> bool:
    if not isinstance(func_string, str) or pd.isna(func_string):
        return False
    return target_function in [f.strip() for f in func_string.split(';')]


def redact_string_descriptions(input_str: str, max_length: int = 30) -> str:
    if not isinstance(input_str, str):
        return str(input_str)
    if len(input_str) <= max_length:
        return input_str
    return input_str[:max_length - 3] + '...'


def _coerce_font_size(value, default: int) -> int:
    """Parse an Appearance Settings font-size input, falling back to `default`
    for blank/invalid/non-positive values (e.g. an emptied number input)."""
    try:
        size = int(float(value))
    except (TypeError, ValueError):
        return default
    return size if size > 0 else default


# ---------------------------------------------------------------------------
# DataAnalysisState – holds all computed intermediate data
# ---------------------------------------------------------------------------

class DataAnalysisState:
    """Holds all intermediate computation results (replaces widget state)."""

    def __init__(self, merged_df, group_data_dict, protein_dict, params: dict):
        self.merged_df = merged_df.copy()
        self.group_data_dict = group_data_dict
        self.protein_dict = protein_dict

        # "Strip protein name" toggle (Appearance Settings). Default on: drop the
        # trailing UniProt descriptors so plot labels use the short name. Every
        # display label (Description, selected_protein_names, protein_reps_dict)
        # derives from protein_dict[pid]['name'], so cleaning it here propagates
        # everywhere. Applied to a copy so the on-disk protein_dict is untouched.
        self.strip_protein_name: bool = params.get('strip_protein_name', True)
        if self.strip_protein_name and self.protein_dict:
            self.protein_dict = {
                pid: {**info, 'name': _clean_protein_name(info.get('name', pid)) or pid}
                for pid, info in self.protein_dict.items()
            }

        # Widget-state equivalents from params
        self.selected_groups: list = params.get('selected_groups', list(group_data_dict.keys()))
        self.selected_proteins_raw: list = params.get('selected_proteins', [NO_PROTEIN_FILTER])
        self.selected_functions_raw: list = params.get('selected_functions', [NO_FUNCTION_FILTER])
        self.plot_type: str = params.get('plot_type', 'Grouped Bar Plots')
        # The Plot Filter dropdown is gone from the UI: the mode is derived from
        # what the user selected. An explicit plot_filter is still honoured so
        # programmatic callers (and the test suite) can pin a mode directly.
        self.plot_filter: str = params.get('plot_filter') or derive_plot_filter(
            self.selected_proteins_raw, self.selected_functions_raw)
        # Which dimensions the user actually restricted. Derived from the raw
        # selections rather than from plot_filter so that an explicitly-passed
        # mode and a derived one behave the same.
        self.protein_filter_active: bool = _selection_is_active(
            self.selected_proteins_raw, NO_PROTEIN_FILTER)
        self.function_filter_active: bool = _selection_is_active(
            self.selected_functions_raw, NO_FUNCTION_FILTER)
        self.abs_or_count: str = params.get('abs_or_count', 'Abundance')
        self.metric_type: str = params.get('metric_type', 'Absolute')
        self.orientation: str = params.get('orientation', 'By Sample')
        self.log_transform: bool = params.get('log_transform', False)
        self.plot_minor: bool = params.get('plot_minor', False)
        # Overlay group-comparison significance markers on bar charts.
        self.show_significance: bool = params.get('show_significance', False)
        # Which test drives the markers: 'tukey' (ANOVA + Tukey HSD, equal-variance
        # lab convention) or 'games-howell' (Welch ANOVA + Games–Howell, robust to
        # the unequal variance typical of abundance data). See services/stats.py.
        self.significance_method: str = params.get('significance_method', 'tukey')
        self.color_scheme: str = params.get('color_scheme', 'HSV')
        self.xlabel: str = params.get('xlabel', '')
        self.ylabel: str = params.get('ylabel', '')
        self.legend_title: str = params.get('legend_title', '')
        self.plot_title: str = params.get('plot_title', '')
        # Appearance Settings font-size overrides (Auto = the historical hardcoded
        # sizes baked into plotter.py). Coerced defensively since these arrive as
        # raw JSON from the client and an empty/non-numeric value should fall back
        # to the default rather than blow up figure generation.
        self.font_size_xlabel: int = _coerce_font_size(params.get('font_size_xlabel'), 18)
        self.font_size_ylabel: int = _coerce_font_size(params.get('font_size_ylabel'), 18)
        self.font_size_legend_title: int = _coerce_font_size(params.get('font_size_legend_title'), 14)
        # Tick labels (the numbers/categories along each axis) are a distinct
        # element from the axis title above, so they get their own size —
        # notably the correlation plots show numeric ticks on both axes.
        self.font_size_xtick: int = _coerce_font_size(params.get('font_size_xtick'), 16)
        self.font_size_ytick: int = _coerce_font_size(params.get('font_size_ytick'), 16)
        # Data-value labels drawn above bars/stacked-bars, and the percent text
        # inside pie slices — one shared size for every "number on the chart"
        # element that isn't an axis, title, or legend.
        self.font_size_value_label: int = _coerce_font_size(params.get('font_size_value_label'), 12)
        self.font_size_plot_title: int = _coerce_font_size(params.get('font_size_plot_title'), 18)
        self.correlation_type: str = params.get('correlation_type', 'Pearson')

        # Derived
        self.avg_columns = [c for c in self.merged_df.columns if c.startswith('Avg_')]
        self.stripped_columns = [c.replace('Avg_', '') for c in self.avg_columns]
        self.use_count = self.abs_or_count == 'Count'
        self.is_relative = 'relative' in self.metric_type.lower()

        if self.use_count:
            self.value_prefix = 'Count_'
            self.rel_prefix = 'Rel_Count_'
            self.metric_name = 'Unique Peptide Count'
            self.num_format = ',.0f'
        else:
            self.value_prefix = 'Avg_'
            self.rel_prefix = 'Rel_Avg_'
            self.metric_name = 'Summed Abundance'
            self.num_format = ',.2e'

        self.value_cols = {g: f'{self.value_prefix}{g}' for g in self.selected_groups}
        self.rel_cols = {g: f'{self.rel_prefix}{g}' for g in self.selected_groups}

        # Non-fatal, user-facing notices collected while building a plot (e.g. a
        # zero/negative value dropped under log transform). Surfaced alongside
        # the figure by generate_plot() rather than raised as an error.
        self.warnings: list = []

        # Computed by pipeline steps
        self.filtered_df: pd.DataFrame | None = None
        self.total_peptide_results_dict: dict = {}
        self.protein_df: pd.DataFrame | None = None
        self.sum_df: pd.DataFrame | None = None
        self.protein_sample_distribution_dict: dict = {}
        self.protein_count_bysample_dict: dict = {}
        self.function_df: pd.DataFrame | None = None
        self.function_distribution_dict: dict = {}
        self.function_abundance_totals_dict: dict = {}
        self.function_count_totals_dict: dict = {}
        self.unique_function_abundance_dict: dict = {}
        self.unique_function_counts_dict: dict = {}
        self.function_group_metrics_dict: dict = {}
        self.function_sem_dict: dict = {}  # {group: {fn: {'abundance_sem':.., 'count_sem':..}}}
        # Per-replicate value vectors for group-comparison statistics (ANOVA/Tukey):
        #   {category_or_group: {group: {'abundance': [...], 'count': [...]}}}
        self.total_reps_dict: dict = {}      # {group: {'abundance':[...], 'count':[...]}}
        self.function_reps_dict: dict = {}   # {fn: {group: {...}}}
        self.protein_reps_dict: dict = {}    # {protein_name: {group: {...}}}
        self.function_unique_peptide_counts: dict = {}  # {fn: total unique peptides for that fn}
        self.abundance_count_by_sample_dict: dict = {}
        self.selected_proteins: list = []
        self.selected_functions: list = []
        self.all_proteins: list = []
        self.all_functions: list = []
        self.function_color_map: dict = {}

    @property
    def selected_protein_names(self) -> list:
        """Return selected protein IDs converted to display names (for protein_df/psdd lookups)."""
        return [self.protein_dict.get(pid, {}).get('name', pid) for pid in self.selected_proteins]

    # ------------------------------------------------------------------
    # Step 1: Resolve selected proteins / functions
    # ------------------------------------------------------------------

    def _resolve_all_proteins(self):
        """Build all_proteins list sorted by abundance."""
        protein_abundance = _compute_protein_abundance(self.merged_df, self.avg_columns)

        self.all_proteins = sorted(
            self.protein_dict.keys(),
            key=lambda p: protein_abundance.get(p, 0),
            reverse=True,
        )

    def _resolve_selected_proteins(self):
        """Determine which proteins are actually selected."""
        self._resolve_all_proteins()
        names_to_ids = {self.protein_dict[pid].get('name', pid): pid for pid in self.protein_dict}

        # "All Proteins", or no protein restriction at all, both resolve to every
        # protein: downstream lookups (protein_df, reps dicts) always need a
        # concrete list — whether the plot breaks out per protein is decided by
        # plot_filter, not here.
        if (ALL_PROTEINS in self.selected_proteins_raw
                or not self.protein_filter_active):
            self.selected_proteins = list(self.all_proteins)
            return

        resolved = []
        for sel in self.selected_proteins_raw:
            if sel == NO_PROTEIN_FILTER:
                continue
            if sel in self.protein_dict:
                resolved.append(sel)
            elif sel in names_to_ids:
                resolved.append(names_to_ids[sel])
            else:
                resolved.append(sel)
        self.selected_proteins = resolved or list(self.all_proteins)

    def _resolve_all_functions(self):
        function_totals = _compute_function_counts(self.merged_df)
        self.all_functions = [
            f for f, _ in sorted(function_totals.items(), key=lambda x: x[1], reverse=True)
        ]

    def _resolve_selected_functions(self):
        self._resolve_all_functions()
        broader = {'Functional Peptides', 'Non-Functional Peptides', 'Other Functions'}

        if self.plot_filter == 'Functional vs Non-Functional Peptides':
            self.selected_functions = ['Functional Peptides', 'Non-Functional Peptides']
        elif (ALL_FUNCTIONAL in self.selected_functions_raw
                or not self.function_filter_active):
            # "All Functional Peptides", or no function restriction: same as for
            # proteins, resolve to the full list for downstream lookups.
            self.selected_functions = [f for f in self.all_functions if f not in broader]
        else:
            self.selected_functions = [
                f for f in self.selected_functions_raw
                if f not in broader and f != NO_FUNCTION_FILTER
            ]

        # Generate function color map
        colors = get_color_sequence(len(self.all_functions), self.color_scheme)
        self.function_color_map = {func: color for func, color in zip(self.all_functions, colors)}
        self.function_color_map.setdefault('Functional Peptides', '#4CAF50')
        self.function_color_map.setdefault('Non-Functional Peptides', '#9E9E9E')
        self.function_color_map.setdefault('Other Functions', '#808080')

    # ------------------------------------------------------------------
    # Step 2: Filter dataframe + total peptide results
    # ------------------------------------------------------------------

    def filter_dataframe(self):
        df = self.merged_df.copy()
        protein_mask = pd.Series(True, index=df.index)
        function_mask = pd.Series(True, index=df.index)

        # Protein filter – match by ID against the 'Protein' column (semicolon-separated IDs).
        # When plot_minor=True, skip the protein restriction so all protein rows remain in
        # filtered_df; process_protein_data() will aggregate non-selected ones into "Other Proteins".
        if (self.plot_filter in ('Selected Protein(s)', 'Both',
                                 'Functional vs Non-Functional Peptides')
                and self.protein_filter_active
                and 'Protein' in df.columns
                and not self.plot_minor):
            if ALL_PROTEINS not in self.selected_proteins_raw:
                sel_ids = set(self.selected_proteins)
                protein_mask = df['Protein'].fillna('').apply(
                    lambda x: bool(sel_ids.intersection(p.strip() for p in x.split(';')))
                )

        # Function filter
        has_function = 'function' in df.columns and not df['function'].isna().all()
        if has_function and self.plot_filter in ('Selected Function(s)', 'Both', 'Functional vs Non-Functional Peptides'):
            if self.plot_filter == 'Functional vs Non-Functional Peptides':
                pass  # don't filter here
            elif self.plot_minor:
                # When grouping minor functions, keep all functional rows so
                # create_function_df() can aggregate non-selected ones into "Other Functions".
                function_mask = df['function'].notna() & (df['function'] != '')
            elif 'All Functional Peptides' in self.selected_functions_raw:
                function_mask = df['function'].notna() & (df['function'] != '')
            elif 'Non-Functional Peptides' in self.selected_functions_raw:
                function_mask = df['function'].isna()
            else:
                sel_funcs = self.selected_functions
                function_mask = df['function'].apply(
                    lambda x: any(contains_function(x, f) for f in sel_funcs)
                )

        if self.plot_filter == 'Both':
            combined = protein_mask & function_mask
        elif self.plot_filter == 'Selected Protein(s)':
            combined = protein_mask
        elif self.plot_filter == 'Functional vs Non-Functional Peptides':
            # function_mask stays all-True here (the split is made downstream, not
            # by filtering); the protein restriction, if any, still applies.
            combined = protein_mask & function_mask
        elif self.plot_filter == 'Selected Function(s)':
            combined = function_mask
        else:
            combined = pd.Series(True, index=df.index)

        self.filtered_df = df[combined]

        # Build total_peptide_results_dict
        total_results = {}
        for group_name in self.selected_groups:
            col = f'Avg_{group_name}'
            if col not in self.filtered_df.columns:
                continue
            temp = self.filtered_df[['Unique Peptide ID', col]].copy()
            temp[col] = pd.to_numeric(temp[col], errors='coerce')
            nonzero = temp[(temp[col].notna()) & (temp[col] != 0)]
            total_abundance = nonzero[col].sum()
            unique_count = nonzero['Unique Peptide ID'].nunique()

            # SEM + per-replicate vectors (need replicate data if available)
            rep_vals, rep_counts = self.replicate_values(self.filtered_df, group_name)
            abundance_sem = (
                float(np.std(rep_vals, ddof=1) / np.sqrt(len(rep_vals))) if len(rep_vals) > 1 else 0.0
            )
            count_sem = (
                float(np.std(rep_counts, ddof=1) / np.sqrt(len(rep_counts))) if len(rep_counts) > 1 else 0.0
            )
            self.total_reps_dict[group_name] = {'abundance': rep_vals, 'count': rep_counts}

            total_results[group_name] = {
                'total_Abundance': float(total_abundance),
                'unique_peptides': int(unique_count),
                'abundance_sem': abundance_sem,
                'count_sem': count_sem,
                'relative_Abundance': 0.0,
                'relative_peptides': 0.0,
            }

        # Calculate relatives
        total_ab = sum(v['total_Abundance'] for v in total_results.values())
        total_pep = sum(v['unique_peptides'] for v in total_results.values())
        for g, v in total_results.items():
            v['relative_Abundance'] = v['total_Abundance'] / total_ab * 100 if total_ab else 0
            v['relative_peptides'] = v['unique_peptides'] / total_pep * 100 if total_pep else 0

        self.total_peptide_results_dict = total_results

        # Also build abundance_count_by_sample_dict (same as total_results for now)
        self.abundance_count_by_sample_dict = {
            g: {
                'total_Abundance': v['total_Abundance'],
                'unique_peptides': v['unique_peptides'],
                'relative_Abundance': v['relative_Abundance'],
                'relative_peptides': v['relative_peptides'],
            }
            for g, v in total_results.items()
        }

    # ------------------------------------------------------------------
    # Replicate-based SEM helper (shared by function/protein bar charts)
    # ------------------------------------------------------------------

    def replicate_values(self, sub_df, group_name):
        """Per-replicate summed abundance and unique-peptide count for the peptides
        in ``sub_df``, one value per replicate column of ``group_name``.

        These per-replicate vectors are the raw material for both the SEM error
        bars and the group-comparison statistics (ANOVA / Tukey HSD): each
        replicate contributes one summed abundance and one presence-count.

        Returns (rep_abundances, rep_counts) as lists. Returns ([], []) when the
        group has fewer than two replicate columns (the no-replicate input form).
        """
        gd = self.group_data_dict.get(group_name)
        if not isinstance(gd, list) or len(gd) < 2 or sub_df is None or len(sub_df) == 0:
            return [], []
        rep_vals, rep_counts = [], []
        has_id = 'Unique Peptide ID' in sub_df.columns
        for rep_col in gd:
            if rep_col not in sub_df.columns:
                continue
            vals = pd.to_numeric(sub_df[rep_col], errors='coerce')
            rep_vals.append(float(vals.fillna(0).sum()))
            nz = vals.notna() & (vals > 0)
            if has_id:
                rep_counts.append(int(sub_df.loc[nz, 'Unique Peptide ID'].nunique()))
            else:
                rep_counts.append(int(nz.sum()))
        return rep_vals, rep_counts

    def replicate_sems(self, sub_df, group_name):
        """Standard error of the mean across a group's replicate columns, for the
        abundance and unique-peptide count of the peptides in ``sub_df``.

        Same replicate-SEM definition used for the sample-totals bar chart (see
        ``filter_dataframe``); applied to a category subset so grouped bars
        oriented By Function / By Protein carry error bars too. SEM =
        std(ddof=1)/sqrt(n). Returns (0.0, 0.0) for the no-replicate input form.
        """
        rep_vals, rep_counts = self.replicate_values(sub_df, group_name)
        ab_sem = float(np.std(rep_vals, ddof=1) / np.sqrt(len(rep_vals))) if len(rep_vals) > 1 else 0.0
        ct_sem = float(np.std(rep_counts, ddof=1) / np.sqrt(len(rep_counts))) if len(rep_counts) > 1 else 0.0
        return ab_sem, ct_sem

    # ------------------------------------------------------------------
    # Step 3: Calculate bioactive function data
    # ------------------------------------------------------------------

    def calculate_bioactive_data(self):
        df = self.filtered_df
        if df is None or 'function' not in df.columns:
            return

        all_fns = (
            ['Functional Peptides', 'Non-Functional Peptides']
            if self.plot_filter == 'Functional vs Non-Functional Peptides'
            else self.all_functions
        )

        unique_fn_ab: dict = {}
        unique_fn_ct: dict = {}
        fn_ab_totals: dict = {}
        fn_ct_totals: dict = {}
        fn_group_metrics: dict = {}

        # Per-function unique peptide counts across all groups (prevents double counting in totals).
        # A peptide with functions "A;B" is counted in BOTH A and B individually,
        # but total unique peptide count uses nunique() on the full set.
        fn_unique_peptide_counts: dict = {}
        # Cache each function's row-membership mask on the full filtered df so we can
        # reuse it below to compute replicate-based SEM (needs the raw replicate
        # columns, which `temp`/`nonzero` do not carry).
        fn_df_masks: dict = {}
        for fn in all_fns:
            if self.plot_filter == 'Functional vs Non-Functional Peptides':
                mask = df['function'].notna() if fn == 'Functional Peptides' else df['function'].isna()
            else:
                mask = df['function'].apply(lambda x: contains_function(x, fn))
            fn_df_masks[fn] = mask
            fn_rows = df[mask]
            fn_unique_peptide_counts[fn] = (
                int(fn_rows['Unique Peptide ID'].nunique())
                if 'Unique Peptide ID' in fn_rows.columns else 0
            )
        self.function_unique_peptide_counts = fn_unique_peptide_counts

        # SEM per (group, function): {group: {fn: {'abundance_sem':.., 'count_sem':..}}}
        fn_sem: dict = {}

        for group_name in self.selected_groups:
            col = f'Avg_{group_name}'
            if col not in df.columns:
                continue
            temp = df[['Unique Peptide ID', 'function', col]].copy()
            temp[col] = pd.to_numeric(temp[col], errors='coerce')
            nonzero = temp[(temp[col] != 0) & temp[col].notna()]
            unique_fn_ab.setdefault(group_name, {})
            unique_fn_ct.setdefault(group_name, {})
            fn_sem.setdefault(group_name, {})

            for fn in all_fns:
                if self.plot_filter == 'Functional vs Non-Functional Peptides':
                    mask = nonzero['function'].notna() if fn == 'Functional Peptides' else nonzero['function'].isna()
                else:
                    mask = nonzero['function'].apply(lambda x: contains_function(x, fn))

                fn_rows = nonzero[mask]
                ab = float(fn_rows[col].sum())
                ct = int(fn_rows['Unique Peptide ID'].nunique())
                unique_fn_ab[group_name][fn] = ab
                unique_fn_ct[group_name][fn] = ct

                # Replicate vectors for this (group, function), from the raw
                # replicate columns of the peptides annotated to this function —
                # drive both the SEM error bar and the group-comparison stats.
                fn_sub = df[fn_df_masks[fn]]
                rep_ab, rep_ct = self.replicate_values(fn_sub, group_name)
                ab_sem = float(np.std(rep_ab, ddof=1) / np.sqrt(len(rep_ab))) if len(rep_ab) > 1 else 0.0
                ct_sem = float(np.std(rep_ct, ddof=1) / np.sqrt(len(rep_ct))) if len(rep_ct) > 1 else 0.0
                fn_sem[group_name][fn] = {'abundance_sem': ab_sem, 'count_sem': ct_sem}
                self.function_reps_dict.setdefault(fn, {})[group_name] = {
                    'abundance': rep_ab, 'count': rep_ct,
                }

                fn_group_metrics.setdefault(fn, {})
                fn_group_metrics[fn][group_name] = {
                    'Abundance': ab,
                    'count': ct,
                    'rel_Abundance': 0.0,
                    'rel_count': 0.0,
                }

            # CRITICAL: Use actual unique functional-peptide counts (not sum of per-function
            # counts) to prevent double-counting peptides with multiple function annotations.
            # E.g. if a peptide has "A;B", summing fn_ct[A]+fn_ct[B] counts it twice.
            if self.plot_filter == 'Functional vs Non-Functional Peptides':
                # Count all non-null rows (functional) and null rows separately
                functional_nonzero = nonzero[nonzero['function'].notna()]
                nonfunctional_nonzero = nonzero[nonzero['function'].isna()]
                fn_ab_totals[group_name] = float(functional_nonzero[col].sum()) + float(nonfunctional_nonzero[col].sum())
                fn_ct_totals[group_name] = int(nonzero['Unique Peptide ID'].nunique())
            else:
                # Rows with any functional annotation (may have multi-function peptides)
                func_mask = nonzero['function'].notna() & (nonzero['function'].astype(str).str.strip() != '')
                functional_nonzero = nonzero[func_mask]
                # Total functional abundance: sum each unique peptide once (not once per function)
                fn_ab_totals[group_name] = float(
                    functional_nonzero.drop_duplicates('Unique Peptide ID')[col].sum()
                ) if 'Unique Peptide ID' in functional_nonzero.columns else float(functional_nonzero[col].sum())
                # Total functional peptide count: unique peptides with any function
                fn_ct_totals[group_name] = int(functional_nonzero['Unique Peptide ID'].nunique()) if 'Unique Peptide ID' in functional_nonzero.columns else 0

        # Relative values: function's contribution to each group's functional total
        for fn, gdict in fn_group_metrics.items():
            for gname, vals in gdict.items():
                tot_ab = fn_ab_totals.get(gname, 0)
                tot_ct = fn_ct_totals.get(gname, 0)
                vals['rel_Abundance'] = vals['Abundance'] / tot_ab * 100 if tot_ab else 0
                vals['rel_count'] = vals['count'] / tot_ct * 100 if tot_ct else 0

        self.unique_function_abundance_dict = unique_fn_ab
        self.unique_function_counts_dict = unique_fn_ct
        self.function_abundance_totals_dict = fn_ab_totals
        self.function_count_totals_dict = fn_ct_totals
        self.function_group_metrics_dict = fn_group_metrics
        self.function_sem_dict = fn_sem

    # ------------------------------------------------------------------
    # Step 4: Build function DataFrame
    # ------------------------------------------------------------------

    def create_function_df(self):
        if self.filtered_df is None:
            return

        if self.plot_filter == 'Functional vs Non-Functional Peptides':
            all_fns = ['Functional Peptides', 'Non-Functional Peptides']
        else:
            broader = {'Functional Peptides', 'Non-Functional Peptides', 'Other Functions'}
            all_fns = [f for f in self.all_functions if f not in broader]

        rows = []
        for fn in all_fns:
            row = {'Description': fn}
            for g in self.selected_groups:
                ab = self.unique_function_abundance_dict.get(g, {}).get(fn, 0)
                ct = self.unique_function_counts_dict.get(g, {}).get(fn, 0)
                sem = self.function_sem_dict.get(g, {}).get(fn, {})
                row[f'Avg_{g}'] = ab
                row[f'Rel_Avg_{g}'] = 0.0
                row[f'Count_{g}'] = ct
                row[f'Rel_Count_{g}'] = 0.0
                row[f'SEM_Avg_{g}'] = float(sem.get('abundance_sem', 0.0))
                row[f'SEM_Count_{g}'] = float(sem.get('count_sem', 0.0))
            rows.append(row)

        if not rows:
            self.function_df = pd.DataFrame()
            return

        self.function_df = pd.DataFrame(rows)

        for g in self.selected_groups:
            tot_ab = self.function_df[f'Avg_{g}'].sum()
            tot_ct = self.function_df[f'Count_{g}'].sum()
            if tot_ab > 0:
                self.function_df[f'Rel_Avg_{g}'] = (self.function_df[f'Avg_{g}'] / tot_ab * 100).round(6)
            if tot_ct > 0:
                self.function_df[f'Rel_Count_{g}'] = (self.function_df[f'Count_{g}'] / tot_ct * 100).round(6)

        ab_cols = [f'Avg_{g}' for g in self.selected_groups]
        self.function_df['avg_abundance_all'] = (
            self.function_df[ab_cols].sum(axis=1) /
            max(self.function_df[ab_cols].sum().sum(), 1e-10) * 100
        ).round(6)
        self.function_df = self.function_df.sort_values('avg_abundance_all', ascending=False)

        # Handle "Other Functions"
        if self.plot_minor and self.plot_filter not in ('Functional vs Non-Functional Peptides',):
            selected_fn_set = set(self.selected_functions)
            minor_rows = self.function_df[~self.function_df['Description'].isin(selected_fn_set)]
            main_rows = self.function_df[self.function_df['Description'].isin(selected_fn_set)]

            if not minor_rows.empty:
                minor_row = {'Description': 'Other Functions'}
                for g in self.selected_groups:
                    minor_row[f'Avg_{g}'] = minor_rows[f'Avg_{g}'].sum()
                    minor_row[f'Rel_Avg_{g}'] = minor_rows[f'Rel_Avg_{g}'].sum()
                    minor_row[f'Count_{g}'] = minor_rows[f'Count_{g}'].sum()
                    minor_row[f'Rel_Count_{g}'] = minor_rows[f'Rel_Count_{g}'].sum()
                    # SEM does not sum across the pooled minor functions; leave the
                    # aggregate "Other Functions" bar without an error bar rather
                    # than report a statistically meaningless summed SEM.
                    minor_row[f'SEM_Avg_{g}'] = 0.0
                    minor_row[f'SEM_Count_{g}'] = 0.0
                minor_row['avg_abundance_all'] = minor_rows['avg_abundance_all'].sum()
                self.function_df = pd.concat(
                    [main_rows, pd.DataFrame([minor_row])],
                    ignore_index=True,
                )
                if 'Other Functions' not in self.selected_functions:
                    self.selected_functions = list(self.selected_functions) + ['Other Functions']

    # ------------------------------------------------------------------
    # Step 5: Build protein DataFrame
    # ------------------------------------------------------------------

    def process_protein_data(self):
        if self.filtered_df is None or not self.protein_dict:
            return

        df = self.filtered_df.copy()
        abundance_cols = [f'Avg_{g}' for g in self.selected_groups]
        available_ab = [c for c in abundance_cols if c in df.columns]

        if not available_ab or 'Protein' not in df.columns:
            return

        # Calculate protein totals – track per-group peptide sets for Count columns.
        #
        # Vectorized: explode() turns "one row per peptide" into "one row per
        # (peptide, protein-id) pair" without a Python iterrows() pass over
        # every peptide row. What's left is a Python loop, but over distinct
        # protein IDs (a groupby) rather than over every row -- on a dataset
        # with many peptides per protein that's an order-of-magnitude fewer
        # iterations, and each iteration's work (sums, boolean masks, set
        # construction) is a vectorized pandas op rather than per-cell access.
        pep_id_col = 'Unique Peptide ID' if 'Unique Peptide ID' in df.columns else None
        pid_parts = df['Protein'].astype(str).str.split(';')
        exploded = df.assign(_pid=pid_parts).explode('_pid')
        exploded['_pid'] = exploded['_pid'].str.strip()
        exploded = exploded[exploded['_pid'] != '']

        protein_data = {}
        if not exploded.empty:
            for pid, g in exploded.groupby('_pid', sort=False):
                # A plain numpy bool array indexes positionally, sidestepping any
                # ambiguity from duplicate index labels (explode() repeats the
                # original row's index for each protein ID it produced).
                pep_ids_here = g[pep_id_col] if pep_id_col else pd.Series([''] * len(g), index=g.index)
                abundance = {}
                peptides_by_group = {}
                for grp in self.selected_groups:
                    col = f'Avg_{grp}'
                    if col in g.columns:
                        col_vals = pd.to_numeric(g[col], errors='coerce')
                        abundance[grp] = float(col_vals.fillna(0).sum())
                        peptides_by_group[grp] = set(pep_ids_here[(col_vals > 0).to_numpy()])
                    else:
                        abundance[grp] = 0.0
                        peptides_by_group[grp] = set()
                protein_data[pid] = {
                    'peptides': set(pep_ids_here),
                    'peptides_by_group': peptides_by_group,
                    'abundance': abundance,
                    'row_idx': list(g.index),  # filtered_df rows mapped to this protein (for SEM)
                }

        rows = []
        for pid, pdata in protein_data.items():
            name = self.protein_dict.get(pid, {}).get('name', pid)
            row = {'Protein': pid, 'Description': name}
            # Rows of filtered_df mapped to this protein — used for replicate SEM/stats.
            prot_rows = df.loc[pdata['row_idx']] if pdata['row_idx'] else None
            for g in self.selected_groups:
                ab = pdata['abundance'][g]
                row[f'Avg_{g}'] = ab
                row[f'Count_{g}'] = len(pdata['peptides_by_group'][g])
                rep_ab, rep_ct = self.replicate_values(prot_rows, g)
                row[f'SEM_Avg_{g}'] = (
                    float(np.std(rep_ab, ddof=1) / np.sqrt(len(rep_ab))) if len(rep_ab) > 1 else 0.0
                )
                row[f'SEM_Count_{g}'] = (
                    float(np.std(rep_ct, ddof=1) / np.sqrt(len(rep_ct))) if len(rep_ct) > 1 else 0.0
                )
                self.protein_reps_dict.setdefault(name, {})[g] = {
                    'abundance': rep_ab, 'count': rep_ct,
                }
            row['unique_peptide_count'] = len(pdata['peptides'])
            rows.append(row)

        if not rows:
            self.protein_df = pd.DataFrame()
            return

        self.protein_df = pd.DataFrame(rows)

        # Relative abundances (protein's contribution to each group total)
        for g in self.selected_groups:
            ab_col = f'Avg_{g}'
            ct_col = f'Count_{g}'
            tot_ab = self.protein_df[ab_col].sum()
            tot_ct = self.protein_df[ct_col].sum()
            self.protein_df[f'Rel_Avg_{g}'] = (self.protein_df[ab_col] / tot_ab * 100).round(6) if tot_ab else 0.0
            self.protein_df[f'Rel_Count_{g}'] = (self.protein_df[ct_col] / tot_ct * 100).round(6) if tot_ct else 0.0

        # sum_df
        self.sum_df = pd.DataFrame({
            'Sample': [f'Avg_{g}' for g in self.selected_groups],
            'Total_Sum': [self.protein_df[f'Avg_{g}'].sum() for g in self.selected_groups],
        })

        # Sort by total abundance
        ab_cols = [f'Avg_{g}' for g in self.selected_groups]
        total_sum = max(self.protein_df[ab_cols].sum().sum(), 1e-10)
        self.protein_df['avg_abundance_all'] = (
            self.protein_df[ab_cols].sum(axis=1) / total_sum * 100
        ).round(2)
        self.protein_df = self.protein_df.sort_values('avg_abundance_all', ascending=False)

        # ── Handle Other Proteins ──────────────────────────────────────────────
        # When plot_minor=True and specific proteins are selected (not "All"),
        # aggregate non-selected proteins into a single "Other Proteins" entry.
        if (self.plot_minor
                and 'All Proteins (No Filter)' not in self.selected_proteins_raw
                and self.selected_proteins):
            selected_names_set = set(self.selected_protein_names)
            minor_rows = self.protein_df[~self.protein_df['Description'].isin(selected_names_set)]
            main_rows = self.protein_df[self.protein_df['Description'].isin(selected_names_set)]

            if not minor_rows.empty:
                minor_entry = {'Protein': 'Other Proteins', 'Description': 'Other Proteins'}
                for g in self.selected_groups:
                    ab_col = f'Avg_{g}'
                    ct_col = f'Count_{g}'
                    minor_entry[ab_col] = float(minor_rows[ab_col].sum()) if ab_col in minor_rows.columns else 0.0
                    minor_entry[ct_col] = float(minor_rows[ct_col].sum()) if ct_col in minor_rows.columns else 0.0
                    # No meaningful summed SEM for the pooled "Other Proteins" bar.
                    minor_entry[f'SEM_Avg_{g}'] = 0.0
                    minor_entry[f'SEM_Count_{g}'] = 0.0
                minor_entry['unique_peptide_count'] = (
                    int(minor_rows['unique_peptide_count'].sum())
                    if 'unique_peptide_count' in minor_rows.columns else 0
                )

                # Keep only selected proteins + "Other Proteins"
                self.protein_df = pd.concat(
                    [main_rows, pd.DataFrame([minor_entry])], ignore_index=True
                )

                # Recompute relative columns after aggregation
                for g in self.selected_groups:
                    ab_col = f'Avg_{g}'
                    ct_col = f'Count_{g}'
                    tot_ab = float(self.protein_df[ab_col].sum()) if ab_col in self.protein_df.columns else 0
                    tot_ct = float(self.protein_df[ct_col].sum()) if ct_col in self.protein_df.columns else 0
                    self.protein_df[f'Rel_Avg_{g}'] = (
                        (self.protein_df[ab_col] / tot_ab * 100).round(6) if tot_ab > 0 else 0.0
                    )
                    self.protein_df[f'Rel_Count_{g}'] = (
                        (self.protein_df[ct_col] / tot_ct * 100).round(6) if tot_ct > 0 else 0.0
                    )

                # Recompute avg_abundance_all and sort
                ab_cols_avail = [f'Avg_{g}' for g in self.selected_groups if f'Avg_{g}' in self.protein_df.columns]
                total_sum = max(self.protein_df[ab_cols_avail].sum().sum(), 1e-10)
                self.protein_df['avg_abundance_all'] = (
                    self.protein_df[ab_cols_avail].sum(axis=1) / total_sum * 100
                ).round(2)
                self.protein_df = self.protein_df.sort_values('avg_abundance_all', ascending=False)

                # Update sum_df
                self.sum_df = pd.DataFrame({
                    'Sample': [f'Avg_{g}' for g in self.selected_groups],
                    'Total_Sum': [
                        float(self.protein_df[f'Avg_{g}'].sum())
                        if f'Avg_{g}' in self.protein_df.columns else 0
                        for g in self.selected_groups
                    ],
                })

                # Add "Other Proteins" to selected_proteins list
                if 'Other Proteins' not in self.selected_proteins:
                    self.selected_proteins = list(self.selected_proteins) + ['Other Proteins']

        # Build protein_sample_distribution_dict.
        # 'abundance_relative' = protein's contribution to each group total (for By Sample hover).
        # 'relative'           = each group's contribution to the protein's total (for By Protein stacked bar).
        psdd = {}
        for _, row in self.protein_df.iterrows():
            pname = row['Description']
            tot_ab = sum(float(row.get(f'Avg_{g}', 0)) for g in self.selected_groups)
            tot_ct = sum(float(row.get(f'Count_{g}', 0)) for g in self.selected_groups)
            psdd[pname] = {
                'Abundance': {g: float(row.get(f'Avg_{g}', 0)) for g in self.selected_groups},
                'counts': {g: float(row.get(f'Count_{g}', 0)) for g in self.selected_groups},
                # protein's share of each group's total (By Sample hover / reference)
                'abundance_relative': {g: float(row.get(f'Rel_Avg_{g}', 0)) for g in self.selected_groups},
                'count_relative': {g: float(row.get(f'Rel_Count_{g}', 0)) for g in self.selected_groups},
                # group's share of this protein's total (for By Protein stacked bar)
                'relative': {
                    g: float(row.get(f'Avg_{g}', 0)) / tot_ab * 100 if tot_ab > 0 else 0.0
                    for g in self.selected_groups
                },
                'count_relative_to_protein': {
                    g: float(row.get(f'Count_{g}', 0)) / tot_ct * 100 if tot_ct > 0 else 0.0
                    for g in self.selected_groups
                },
                'values': {g: float(row.get(f'Avg_{g}', 0)) for g in self.selected_groups},
                'unique_peptide_count': int(row.get('unique_peptide_count', 0)),
                'total_Abundance': tot_ab,
                'total_count': tot_ct,
            }
            # Other Proteins gets grey colour in plots
            if pname == 'Other Proteins':
                psdd[pname]['color'] = '#808080'
        self.protein_sample_distribution_dict = psdd

        # Unique peptide counts per group for selected proteins
        self.protein_count_bysample_dict = {
            g: int(self.protein_df[f'Count_{g}'].sum()) if f'Count_{g}' in self.protein_df.columns else 0
            for g in self.selected_groups
        }

    # ------------------------------------------------------------------
    # Step 6: Reorganize by function
    # ------------------------------------------------------------------

    def reorganize_by_function(self):
        if self.function_df is None or self.function_df.empty:
            return

        fdd = {}
        include_fns = set(self.selected_functions)
        if self.plot_minor:
            include_fns.add('Other Functions')

        for _, row in self.function_df.iterrows():
            fn = row['Description']
            if fn not in include_fns:
                continue

            tot_ab = sum(float(row.get(f'Avg_{g}', 0)) for g in self.selected_groups)
            tot_ct = sum(float(row.get(f'Count_{g}', 0)) for g in self.selected_groups)

            fdd[fn] = {
                'Abundance': {g: float(row.get(f'Avg_{g}', 0)) for g in self.selected_groups},
                'counts': {g: float(row.get(f'Count_{g}', 0)) for g in self.selected_groups},
                # replicate-based SEM per group (for grouped-bar error bars)
                'abundance_sem': {g: float(row.get(f'SEM_Avg_{g}', 0)) for g in self.selected_groups},
                'count_sem': {g: float(row.get(f'SEM_Count_{g}', 0)) for g in self.selected_groups},
                # function's share of each group's functional total (By Sample hover/reference)
                'abundance_relative': {g: float(row.get(f'Rel_Avg_{g}', 0)) for g in self.selected_groups},
                'count_relative': {g: float(row.get(f'Rel_Count_{g}', 0)) for g in self.selected_groups},
                # group's share of this function's total (for By Function stacked bar)
                'relative': {
                    g: float(row.get(f'Avg_{g}', 0)) / tot_ab * 100 if tot_ab > 0 else 0.0
                    for g in self.selected_groups
                },
                'count_relative_to_function': {
                    g: float(row.get(f'Count_{g}', 0)) / tot_ct * 100 if tot_ct > 0 else 0.0
                    for g in self.selected_groups
                },
                'total_Abundance': tot_ab,
                'total_count': tot_ct,
                # Unique peptide count for this function (across all groups, not double-counted)
                'unique_peptide_count': self.function_unique_peptide_counts.get(fn, 0),
                'values': {},
            }
            # 'values' alias matches notebook convention
            if self.use_count:
                fdd[fn]['values'] = fdd[fn]['counts']
            else:
                fdd[fn]['values'] = fdd[fn]['Abundance']

        self.function_distribution_dict = fdd

    # ------------------------------------------------------------------
    # Run full pipeline
    # ------------------------------------------------------------------

    def run_pipeline(self):
        self._resolve_selected_proteins()
        self._resolve_selected_functions()
        self.filter_dataframe()
        self.process_protein_data()
        if 'function' in self.filtered_df.columns and not self.filtered_df['function'].isna().all():
            self.calculate_bioactive_data()
            self.create_function_df()
            self.reorganize_by_function()
