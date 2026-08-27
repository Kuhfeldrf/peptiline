"""
Data processing service for the Heatmap Visualization web app.
Extracted from heatmap_visualization.ipynb DataTransformation and HeatmapDataHandler classes.
"""
import io
import re
import sys
import os
import traceback

import numpy as np
import pandas as pd
from django.conf import settings

# heatmap_renderer, fasta_utils, and uniprot_client now live alongside this
# file in heatmap_viz/services/ — no sys.path manipulation needed.

# Raster resolution for rendered heatmaps. 300 dpi is the usual publication
# floor; the preview and the downloaded PNG are the same encoded image.
EXPORT_DPI = 300


def _coerce_font_size(value, default: int) -> int:
    """Parse an Appearance Settings font-size input, falling back to `default`
    for blank/invalid/non-positive values (e.g. an emptied number input)."""
    try:
        size = int(float(value))
    except (TypeError, ValueError):
        return default
    return size if size > 0 else default

# Protein-name display cleaning. Two kinds of clutter are removed:
#   * Trailing UniProt FASTA metadata — "Beta-casein OS=Bos taurus GN=CSN2 …".
#   * A leading UniProt entry-name token — "LACB_BOVIN Beta-lactoglobulin",
#     "B4GT1_BOVIN Beta-1,4-galactosyltransferase 1" (the "CAS_bovine" leader).
# Entry names are always upper-case ID_SPECIES, so the leader regex is upper-case
# only and requires a descriptive name after it — a bare "LACB_BOVIN" (no space)
# and lower-cased names are left untouched, so ordinary names never get clipped.
# Kept in sync verbatim with data_analysis/services/data_processor so the "Strip
# protein name" toggle behaves identically in both apps.
_FASTA_META_RE = re.compile(r'\s+(?:OS|OX|GN|PE|SV)=')
_ENTRY_NAME_LEADER_RE = re.compile(r'^[A-Z0-9]+_[A-Z0-9]+\s+(?=\S)')


def _clean_protein_name(name) -> str:
    """Strip trailing UniProt FASTA metadata and a leading entry-name token."""
    s = _FASTA_META_RE.split(str(name), 1)[0].strip()
    return _ENTRY_NAME_LEADER_RE.sub('', s, count=1).strip()


def _parse_grouped_replicates(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Parse the ``"<col> 'Grouped: (GrpA; GrpB)'"`` replicate-column convention.

    Mirrors ``data_analysis.services.data_processor.process_group_data``: the
    Data-Transformation export tags every replicate abundance column that feeds
    an ``Avg_<var>`` mean with the group(s) it belongs to (see
    ``data_transformation/services/data_combiner.py:update_column_names_with_groups``).
    We rename those columns back to their base name and return a
    ``{grouping_variable: [base_col, ...]}`` map, so the differential-comparison
    track can reach the per-replicate values a pooled SD requires.

    Returns ``(df_renamed, replicates_by_group)``. A single-average export with
    no ``'Grouped:'`` columns yields an empty map and the frame unchanged —
    comparison stays disabled and single-series heatmaps are unaffected.
    """
    grouped_cols = [c for c in df.columns if " 'Grouped:" in str(c)]
    if not grouped_cols:
        return df, {}

    replicates_by_group: dict = {}
    rename_map: dict = {}
    for col in grouped_cols:
        base = col.split(" 'Grouped:")[0].strip()
        match = re.search(r"\((.*?)\)", col)
        if not match:
            continue
        for grp in (g.strip() for g in match.group(1).split(";")):
            if grp:
                replicates_by_group.setdefault(grp, []).append(base)
        rename_map[col] = base

    return df.rename(columns=rename_map), replicates_by_group


# ---------------------------------------------------------------------------
# File loading helpers
# ---------------------------------------------------------------------------

def load_merged_file(file_obj, filename: str) -> tuple:
    """
    Load and basic-validate a merged data CSV/TSV/XLSX file.
    Returns (df, group_data_dict, protein_dict, col_order, error_msg)
    """
    name_lower = filename.lower()
    try:
        content = file_obj.read()
        if name_lower.endswith('.xlsx'):
            df = pd.read_excel(io.BytesIO(content))
        elif name_lower.endswith('.tsv') or name_lower.endswith('.txt'):
            df = pd.read_csv(io.BytesIO(content), sep='\t', low_memory=False)
        else:
            df = pd.read_csv(io.BytesIO(content), low_memory=False)
    except Exception as exc:
        return None, {}, {}, [], str(exc)

    # Excel-sourced exports can carry thousands of fully-blank trailing rows
    # (all-comma/all-empty lines past the real data). Their columns parse as
    # NaN/float while the real rows parse as str, which is what triggers
    # pandas' "mixed types" DtypeWarning — and every downstream iterrows()/
    # apply() pass over the dataframe then wastes time on rows with no data.
    # Same fix used by Data Transformation's loader and Data Analysis's load_file.
    df = df.dropna(how='all')
    df = df[~df.astype(str).apply(lambda row: row.str.strip().eq('').all(), axis=1)]

    df.columns = df.columns.str.strip()

    # Standardize columns
    df, err = _validate_and_standardize_columns(df)
    if err:
        return None, {}, {}, [], err

    # Required columns
    if 'Protein' not in df.columns:
        return None, {}, {}, [], "Missing required column 'Protein'."

    # Parse the replicate-column framework (renames "'Grouped: (…)'" columns to
    # their base names) so each group can carry the per-replicate columns the
    # differential-comparison track needs for a pooled SD. Absent it the map is
    # empty and comparison is simply unavailable.
    df, replicates_by_group = _parse_grouped_replicates(df)

    # Build group_data_dict from Avg_ columns
    avg_columns = [col for col in df.columns if col.startswith('Avg_')]
    if not avg_columns:
        return None, {}, {}, [], "No abundance columns (Avg_*) found in file."

    col_order = [col.replace('Avg_', '') for col in avg_columns]
    group_data_dict = {}
    for i, col in enumerate(avg_columns, 1):
        group_name = col.replace('Avg_', '')
        group_data_dict[str(i)] = {
            'grouping_variable': group_name,
            'abundance_columns': [col],
            'replicate_columns': replicates_by_group.get(group_name, []),
        }

    # Build protein_dict from data
    protein_dict = _build_protein_dict_from_df(df)

    return df, group_data_dict, protein_dict, col_order, None


def _validate_and_standardize_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, str | None]:
    """Standardize start/end/Protein columns. Returns (df, error_str_or_None)."""
    START_STOP_MAP = {
        'Start position': 'start', 'End position': 'end',
        'Peptide start': 'start', 'Peptide end': 'end',
        'Start': 'start', 'End': 'end',
        'protein_start': 'start', 'protein_end': 'end',
        'start': 'start', 'end': 'end',
        'StartPosition': 'start', 'EndPosition': 'end',
        'Peptide Start': 'start', 'Peptide End': 'end',
        'Peptide start position': 'start', 'Peptide end position': 'end',
    }
    PROTEIN_ID_COLUMNS = [
        'Protein', 'Leading proteins', 'Protein Name', 'Protein Accession',
        'Accession Number', 'ProteinGroupId', 'Protein ID', 'Accession',
        'protein_accession',
    ]

    # Case-insensitive matching (headers are already whitespace-stripped by
    # load_merged_file() before this runs). The 'Grouped:'/'Avg_*' markers
    # handled elsewhere in this module are generated internally by the app in
    # a fixed case and are intentionally left out of this normalization.
    start_stop_map_norm = {k.strip().lower(): v for k, v in START_STOP_MAP.items()}
    protein_id_columns_norm = [c.strip().lower() for c in PROTEIN_ID_COLUMNS]

    if 'start' not in df.columns or 'end' not in df.columns:
        for col in list(df.columns):
            norm = col.strip().lower()
            if norm in start_stop_map_norm:
                new_name = start_stop_map_norm[norm]
                if new_name not in df.columns:
                    df = df.rename(columns={col: new_name})

    if 'start' not in df.columns or 'end' not in df.columns:
        return df, (
            "Missing required peptide position data ('start'/'end', or an "
            "equivalent column such as 'Positions in Proteins', 'Peptide start'/"
            "'Peptide end'). Sequence heatmaps require a start and end position "
            "for every peptide and cannot be generated without this data."
        )

    if 'Protein' not in df.columns:
        cols_norm = {c.strip().lower(): c for c in df.columns}
        for norm_name in protein_id_columns_norm:
            if norm_name in cols_norm:
                df = df.rename(columns={cols_norm[norm_name]: 'Protein'})
                break

    if 'Protein' not in df.columns:
        return df, "Missing required column 'Protein'."

    # Normalise protein IDs (strip FASTA pipe notation)
    df = df.copy()
    df['Protein'] = df['Protein'].astype(str)

    def extract_uniprot(x: str) -> str:
        if '|' in x:
            m = re.search(r'\|([A-Z0-9]+)\|', x)
            if m:
                return m.group(1)
        return x

    df['Protein'] = df['Protein'].apply(extract_uniprot)

    # Strip leading/trailing whitespace from string columns that are used as
    # column name components in calculate_abundance() – prevents mismatches.
    for col in ('Unique Peptide ID', 'function'):
        if col in df.columns and df[col].dtype == object:
            df[col] = df[col].str.strip()

    return df, None


def _build_protein_dict_from_df(df: pd.DataFrame) -> dict:
    """Build {protein_id: {name, species, sequence}} from DataFrame columns."""
    protein_dict = {}
    if 'Protein' not in df.columns:
        return protein_dict

    group_cols = ['protein_name', 'protein_species']
    has_info = all(c in df.columns for c in group_cols)

    for protein_id, group in df.groupby('Protein'):
        if not protein_id or pd.isna(protein_id):
            continue
        name = group['protein_name'].iloc[0] if has_info else str(protein_id)
        species = group['protein_species'].iloc[0] if has_info else 'Unknown'
        raw_name = str(name) if pd.notna(name) else str(protein_id)
        protein_dict[protein_id] = {
            'name': _clean_protein_name(name) or str(protein_id) if pd.notna(name) else str(protein_id),
            # Full, unstripped name — surfaced when the user unchecks "Strip
            # protein name" in Appearance Settings.
            'name_raw': raw_name,
            'species': str(species) if pd.notna(species) else 'Unknown',
            'sequence': '',
        }
    return protein_dict


# ---------------------------------------------------------------------------
# FASTA loading
# ---------------------------------------------------------------------------

def load_fasta_file(file_obj, merged_protein_ids: set | None = None) -> tuple[dict, str | None]:
    """
    Parse an uploaded FASTA file and return (protein_dict, error_msg).
    If merged_protein_ids is provided, only include proteins that match.
    """
    try:
        from .fasta_utils import validate_fasta_format, parse_fasta
        from utils.uniprot_client import UniProtClient
        from django.conf import settings as _dj_settings
        _SPEC_TRANSLATE_LIST = getattr(_dj_settings, 'SPEC_TRANSLATE_LIST', [])
    except ImportError:
        pass

    content = file_obj.read()
    if isinstance(content, bytes):
        content = content.decode('utf-8', errors='ignore')

    # Simple FASTA parser fallback
    protein_dict = _parse_fasta_simple(content)

    if merged_protein_ids:
        protein_dict = {k: v for k, v in protein_dict.items() if k in merged_protein_ids}

    return protein_dict, None


def _parse_fasta_simple(content: str) -> dict:
    """Simple FASTA parser returning {protein_id: {name, species, sequence}}."""
    proteins = {}
    current_id = None
    current_name = ''
    current_seq = []

    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_id and current_seq:
                proteins[current_id] = {
                    'name': current_name or current_id,
                    'species': '',
                    'sequence': ''.join(current_seq),
                }
            header = line[1:]
            # Try UniProt format: >sp|P12345|NAME_HUMAN ...
            m = re.match(r'(?:sp|tr)\|([A-Z0-9]+)\|(\S+)\s*(.*)', header)
            if m:
                current_id = m.group(1)
                current_name = _clean_protein_name(m.group(3)) or m.group(2)
            else:
                parts = header.split(None, 1)
                current_id = parts[0]
                current_name = _clean_protein_name(parts[1]) if len(parts) > 1 else current_id
            current_seq = []
        else:
            current_seq.append(line)

    if current_id and current_seq:
        proteins[current_id] = {
            'name': current_name or current_id,
            'species': '',
            'sequence': ''.join(current_seq),
        }
    return proteins


# ---------------------------------------------------------------------------
# UniProt sequence fetching
# ---------------------------------------------------------------------------

def _resolve_mnemonic_accession(protein_id: str, species: str | None = None) -> str | None:
    """Resolve a UniProt *mnemonic* (entry-name prefix, e.g. 'APOA1' or the
    pepnets-style short codes like 'A1AT') to an accession, when `protein_id`
    isn't itself a valid accession.

    UniProt entry names are ``{MNEMONIC}_{ORGANISM_CODE}`` (e.g.
    'APOA1_HUMAN'); the accession-only fetch in `fetch_sequence_from_uniprot`
    can't resolve a bare mnemonic. This searches the `id` field with a
    wildcard organism suffix, restricted to reviewed (Swiss-Prot) entries,
    and ranks matches by organism so results are deterministic when a
    mnemonic exists in multiple species (as most do).
    """
    try:
        import requests
        resp = requests.get(
            'https://rest.uniprot.org/uniprotkb/search',
            params={
                'query': f'id:{protein_id}_* AND reviewed:true',
                'fields': 'accession,id,organism_name',
                'format': 'json',
                'size': 25,
            },
            timeout=10,
        )
        if resp.status_code != 200:
            return None
        results = resp.json().get('results') or []
        if not results:
            return None

        # Rank: species hint (if resolvable to a known UniProt organism code)
        # first, else Human, else the first (UniProt's own relevance order).
        preferred_codes = []
        if species:
            for spec_group in settings.SPEC_TRANSLATE_LIST:
                if species.lower() in [t.lower() for t in spec_group]:
                    preferred_codes.append(spec_group[-1].upper())
                    break
        preferred_codes.append('HUMAN')

        for code in preferred_codes:
            for entry in results:
                if entry.get('uniProtkbId', '').upper().endswith(f'_{code}'):
                    return entry.get('primaryAccession')

        return results[0].get('primaryAccession')
    except Exception:
        return None


def fetch_sequence_from_uniprot(
    protein_id: str, species: str | None = None,
) -> tuple[str | None, int | None]:
    """Fetch protein sequence (+ signal-peptide end, if annotated) from UniProt.

    `protein_id` is normally a UniProt accession (e.g. 'P02647'). If the
    direct accession fetch fails, falls back to resolving `protein_id` as a
    UniProt mnemonic/entry-name prefix (e.g. 'APOA1', or the short codes used
    by some external peptidomics datasets like pepnets) via
    `_resolve_mnemonic_accession`, then retries with the resolved accession.
    `species` (a common name from SPEC_TRANSLATE_LIST, e.g. 'Human') is an
    optional hint used only to disambiguate a mnemonic shared across species.

    Returns (sequence, signal_end); either may be None on failure. signal_end
    is the 1-indexed last residue of the annotated signal peptide, used as the
    default "Strip start sequence" length in Heatmap's Appearance Settings.
    """
    try:
        from utils.uniprot_client import UniProtClient
        client = UniProtClient()
        result = client.fetch_protein_info_with_sequence(protein_id)
        if result:
            _, _, seq, signal_end = result
            if seq:
                return seq, signal_end
    except Exception:
        pass
    try:
        import urllib.request
        url = f'https://www.uniprot.org/uniprot/{protein_id}.fasta'
        with urllib.request.urlopen(url, timeout=10) as resp:
            lines = resp.read().decode('utf-8').splitlines()
            seq = ''.join(l for l in lines if not l.startswith('>'))
            if seq:
                return seq, None
    except Exception:
        pass

    # protein_id wasn't a fetchable accession — try it as a bare mnemonic.
    accession = _resolve_mnemonic_accession(protein_id, species)
    if not accession or accession == protein_id:
        return None, None
    try:
        from utils.uniprot_client import UniProtClient
        client = UniProtClient()
        result = client.fetch_protein_info_with_sequence(accession)
        if result:
            _, _, seq, signal_end = result
            if seq:
                return seq, signal_end
    except Exception:
        pass
    return None, None


# ---------------------------------------------------------------------------
# "Strip start sequence" — UniProt signal-peptide suggestion for the UI
# ---------------------------------------------------------------------------

def get_signal_peptide_suggestions(protein_dict: dict, selected_proteins: list) -> list:
    """
    Return the UniProt-annotated signal-peptide suggestion for each of
    ``selected_proteins``, for the "Strip start sequence" UI to display next
    to the manual residue-count override.

    UniProt's API returns the signal peptide's position (start/end) but not
    its amino-acid letters directly — those are only recoverable by slicing
    the full sequence at that position, which is what this does: the letters
    shown to the user (``signal_sequence``) are ``sequence[:signal_end]``.

    Each entry: {protein_id, name, signal_end, signal_sequence, has_signal}.
    ``signal_end``/``signal_sequence`` are None when UniProt has no "Signal"
    feature annotated for that protein, or its sequence hasn't been fetched
    yet; ``has_signal`` is False in both cases.
    """
    out = []
    for pid in selected_proteins:
        pinfo = protein_dict.get(pid, {})
        signal_end = pinfo.get('signal_end')
        sequence = pinfo.get('sequence', '') or ''
        signal_sequence = None
        if signal_end and sequence:
            try:
                signal_sequence = sequence[:int(signal_end)]
            except (TypeError, ValueError):
                signal_sequence = None
        out.append({
            'protein_id': pid,
            'name': pinfo.get('name', pid),
            'signal_end': signal_end,
            'signal_sequence': signal_sequence,
            'has_signal': signal_end is not None,
        })
    return out


# ---------------------------------------------------------------------------
# Get selector options for the UI
# ---------------------------------------------------------------------------

def get_selector_options(merged_df: pd.DataFrame, group_data_dict: dict, protein_dict: dict, col_order: list) -> dict:
    """Return options for the frontend selectors."""
    # Proteins sorted by occurrence in data
    protein_counts = merged_df['Protein'].value_counts() if 'Protein' in merged_df.columns else pd.Series()
    sorted_proteins = list(protein_counts.index)

    protein_options = []
    for pid in sorted_proteins:
        name = protein_dict.get(pid, {}).get('name', pid)
        has_seq = bool(protein_dict.get(pid, {}).get('sequence', ''))
        protein_options.append({'id': pid, 'label': f"{pid} – {name}", 'has_sequence': has_seq})

    # Variable keys (grouping variables)
    var_key_options = col_order if col_order else [
        v['grouping_variable'] for v in group_data_dict.values()
    ]

    # Check if function data exists
    has_functions = 'function' in merged_df.columns and not merged_df['function'].isna().all()

    # Which grouping variables carry replicate-level ('Grouped:') columns — the
    # differential-comparison track (SMD / log2FC) needs ≥1 replicate per group to
    # form a pooled SD. Surfaced so the UI can populate the Series A/B pickers and
    # show the "comparison available / disabled" banner (Steps 2 & 6).
    var_replicates = {
        v['grouping_variable']: bool(v.get('replicate_columns'))
        for v in group_data_dict.values()
    }
    has_replicates = any(var_replicates.values())

    return {
        'proteins': protein_options,
        'var_keys': var_key_options,
        'has_functions': has_functions,
        'var_replicates': var_replicates,
        'has_replicates': has_replicates,
    }


def get_specific_options(
    merged_df: pd.DataFrame,
    bio_or_pep: str,
    selected_proteins: list = None,
    selected_var_keys: list = None,
) -> list:
    """
    Return options as list of {label, value} dicts for the specific-selector widget.

    bio_or_pep '1' → peptide intervals.
        value = label = "start-end" or "start-end UID"  (matches calculate_abundance() column names)
        Only peptides with at least one non-zero Avg_ value in selected_var_keys are returned.
    bio_or_pep '2' → individual bioactive function names.
        Function column values may be semicolon-delimited; we split them so that
        each individual function appears as a separate selectable option, matching
        the way filter_data_by_selection() does its matching.
    """
    df = merged_df.copy()
    if selected_proteins:
        df = df[df['Protein'].isin(selected_proteins)]

    # For peptide intervals, restrict to rows that have non-zero abundance in at least
    # one of the selected variables — so the dropdown only shows "active" peptides.
    if bio_or_pep == '1' and selected_var_keys:
        avg_cols = [f'Avg_{v}' for v in selected_var_keys if f'Avg_{v}' in df.columns]
        if avg_cols:
            mask = df[avg_cols].apply(pd.to_numeric, errors='coerce').gt(0).any(axis=1)
            df = df[mask]

    # ── Bioactive Functions ────────────────────────────────────────────────────
    if bio_or_pep == '2':
        if 'function' not in df.columns:
            return []
        funcs: set[str] = set()
        for val in df['function'].dropna():
            for part in str(val).split(';'):
                part = part.strip()
                if part and part.lower() != 'nan':
                    funcs.add(part)
        return [{'label': f, 'value': f} for f in sorted(funcs)]

    # ── Peptide Intervals ──────────────────────────────────────────────────────
    if bio_or_pep == '1':
        if 'start' not in df.columns or 'end' not in df.columns:
            return []

        seen: set[str] = set()
        options: list[dict] = []

        for _, row in df.iterrows():
            try:
                start = int(row['start'])
                end = int(row['end'])
            except (ValueError, TypeError):
                continue

            interval = f"{start}-{end}"
            value = interval

            # Append Unique Peptide ID — mirrors calculate_abundance() column naming exactly
            if 'Unique Peptide ID' in df.columns:
                uid_raw = row.get('Unique Peptide ID')
                if pd.notna(uid_raw):
                    uid = str(uid_raw).strip()
                    if uid:
                        value = f"{interval} {uid}"

            if value in seen:
                continue
            seen.add(value)

            options.append({'label': value, 'value': value})

        # Sort by interval start position
        def _sort_key(opt: dict) -> int:
            try:
                return int(opt['value'].split('-')[0])
            except (ValueError, IndexError):
                return 0

        options.sort(key=_sort_key)
        return options

    return []


# ---------------------------------------------------------------------------
# Build available_data_variables_dict
# ---------------------------------------------------------------------------

def build_available_data_variables(
    merged_df: pd.DataFrame,
    protein_dict: dict,
    group_data_dict: dict,
    selected_proteins: list,
    selected_var_keys: list,
    strip_start_sequence: bool = False,
    strip_start_manual: int | None = None,
) -> tuple[dict, list]:
    """
    Build the available_data_variables_dict for the heatmap renderer.

    strip_start_sequence / strip_start_manual implement the "Strip start
    sequence" Appearance-Settings toggle: N residues are physically removed
    from the FRONT of each protein's sequence before rendering, and the
    reduced/mature sequence is what gets plotted. The dataset's peptide
    start/end positions are assumed to already be numbered from that same
    mature protein (position 1 = the first residue after the signal
    peptide), so once the sequence is trimmed no further shift is applied —
    dataset position 1 indexes straight to residue 1 of the reduced
    sequence.

    N, per selected protein: strip_start_manual — a plain residue COUNT
    applied to every selected protein — when given; otherwise falls back to
    each protein's own UniProt-annotated signal-peptide length
    (protein_dict[pid]['signal_end']), else 0 (no trim). A manual value is
    cross-checked against protein_dict[pid]['signal_end'] when that's known
    (see below); for a custom/FASTA-only protein with no UniProt record at
    all, it's simply unverifiable, which is expected and not flagged as an
    error.

    Whether peptide positions also have to MOVE is decided per protein by
    heatmap_renderer.detect_dataset_numbering(), so one checkbox serves both
    real-world cases: a dataset numbered from the mature protein (trim only,
    positions stay put — trimming is what fixes the misalignment) and a
    dataset already aligned to the full precursor (trim AND renumber every
    peptide down by N, so a purely cosmetic strip cannot introduce a shift).

    Returns (available_data_variables_dict, messages).
    """
    # Import from heatmap_renderer (heatmap_viz/services/heatmap_renderer.py)
    try:
        from .heatmap_renderer import (
            export_heatmap_data_to_dict, validate_peptide_alignment, detect_dataset_numbering,
        )
    except ImportError:
        return {}, ['Could not import heatmap_renderer. Check notebook utils path.']

    # Build grouping_var → group_info lookup
    gvar_to_info = {v['grouping_variable']: (k, v) for k, v in group_data_dict.items()}

    available = {}
    messages = []

    for protein_id in selected_proteins:
        # Get protein info
        pinfo = protein_dict.get(protein_id, {})
        protein_name = pinfo.get('name', protein_id)
        protein_name_raw = pinfo.get('name_raw', protein_name)
        protein_species = pinfo.get('species', 'Unknown')
        protein_sequence = pinfo.get('sequence', '')

        # If no sequence, try UniProt
        if not protein_sequence:
            messages.append(f"No sequence found for {protein_id} in protein dictionary or FASTA file.")
            continue

        # Filter merged_df for this protein (needed before the strip block so
        # the numbering detection below can inspect this protein's peptides).
        protein_df = merged_df[merged_df['Protein'] == protein_id].copy()
        if protein_df.empty:
            messages.append(f"No data found for protein {protein_id} in merged data.")
            continue

        # "Strip start sequence" — trim the signal peptide (or a manual
        # residue count) off the FRONT of protein_sequence, so the plotted
        # sequence is the mature protein.
        # A manual override applies to every selected protein; otherwise each
        # protein falls back to its own UniProt-annotated signal-peptide length.
        strip_len = 0
        # UniProt's own signal-peptide annotation for THIS protein (or None if
        # it has none, or was never fetched — e.g. a custom protein list
        # loaded straight from a FASTA with no UniProt lookup at all). Used
        # both for the auto-strip fallback below and to sanity-check a
        # manual override against it, when it's actually known.
        signal_end = pinfo.get('signal_end')
        if strip_start_sequence:
            if strip_start_manual is not None:
                try:
                    strip_len = max(0, int(strip_start_manual))
                except (TypeError, ValueError):
                    strip_len = 0
                # A manual override is one number applied to EVERY selected
                # protein, so it silently goes stale if the user leaves the
                # checkbox on after adding/switching to a protein it wasn't
                # tuned for. Only flag it when there's an actual UniProt
                # value to contradict — a custom/FASTA-only protein simply
                # has nothing to cross-check against, which is expected and
                # not an error.
                if strip_len:
                    if signal_end is None:
                        messages.append(
                            f"{protein_id} has no UniProt-annotated signal peptide on file, so "
                            f"the manual Strip Start Sequence value ({strip_len} residue(s)) "
                            "can't be automatically cross-checked — expected for a custom/"
                            "FASTA-only protein."
                        )
                    elif int(signal_end) != strip_len:
                        messages.append(
                            f"Error: the manual Strip Start Sequence value ({strip_len} "
                            f"residue(s)) does not match {protein_id}'s UniProt-annotated signal "
                            f"peptide length ({signal_end} residue(s)) — verify this is "
                            "intentional before trusting the plotted positions."
                        )
            else:
                if signal_end:
                    try:
                        strip_len = max(0, int(signal_end))
                    except (TypeError, ValueError):
                        strip_len = 0
            strip_len = min(strip_len, len(protein_sequence))

        # Decide whether stripping also has to RENUMBER the peptides. Two
        # different real-world datasets both reach this point:
        #   'mature'    — positions are already numbered from the mature
        #                 protein, so the FASTA's signal peptide is what made
        #                 the plot misaligned; trimming alone fixes it and the
        #                 positions must stay exactly where they are.
        #   'precursor' — positions already line up with the untrimmed FASTA
        #                 (the plot was correct before stripping); the strip is
        #                 purely cosmetic, so every peptide must move down by
        #                 strip_len to stay on the same residue.
        position_offset = 0
        if strip_len:
            numbering, reason = detect_dataset_numbering(
                pinfo.get('sequence', ''), protein_df, strip_len)
            if numbering == 'precursor':
                position_offset = strip_len
                messages.append(
                    f"{protein_id}: peptide positions are numbered against the full precursor "
                    f"({reason}), so they were shifted down by {strip_len} to stay aligned with "
                    "the trimmed sequence."
                )
            else:
                messages.append(
                    f"{protein_id}: peptide positions are numbered from the mature protein "
                    f"({reason}), so {strip_len} residue(s) were trimmed from the sequence "
                    "without moving any positions."
                )
            protein_sequence = protein_sequence[strip_len:]

        is_all_null = (
            'function' not in protein_df.columns
            or protein_df['function'].isna().all()
        )

        # Sanity-check that peptide positions actually land inside the
        # (already-trimmed, if applicable) protein sequence and, when the
        # dataset carries the peptide's own sequence text, that the residues
        # match. Runs once per protein (not per var_key) so a systematic
        # mismatch isn't repeated.
        messages.extend(validate_peptide_alignment(
            protein_sequence, protein_df, position_offset=position_offset,
            protein_id=protein_id,
        ))

        for var_key in selected_var_keys:
            if var_key not in gvar_to_info:
                messages.append(f"Variable key '{var_key}' not found in group data.")
                continue

            group_key, group_info = gvar_to_info[var_key]

            try:
                heatmap_data = export_heatmap_data_to_dict(
                    protein_id, group_key, group_info,
                    protein_sequence, protein_species, protein_name,
                    protein_df, is_all_null, position_offset=position_offset,
                )
            except Exception as exc:
                messages.append(f"Error processing {protein_id}/{var_key}: {exc}")
                continue

            combo_key = f"{protein_id}_{var_key}"
            available[combo_key] = {
                'protein_id': protein_id,
                'protein_sequence': protein_sequence,
                'protein_name': protein_name,
                'protein_name_raw': protein_name_raw,
                'protein_species': protein_species,
                # Per-replicate abundance columns for this group, carried so the
                # differential-comparison track can compute a pooled SD. Empty
                # when the file has no replicate-level ('Grouped:') columns.
                'replicate_columns': group_info.get('replicate_columns', []),
                # Raw peptide rows (with replicate columns) for the differential
                # track; None-safe downstream when comparison is off.
                'peptide_df': heatmap_data.get('peptide_df'),
                # Residues the differential track must also shift peptide_df's
                # (un-renumbered) start/end by — see detect_dataset_numbering.
                'position_offset': position_offset,
                'heatmap_df': heatmap_data.get('heatmap_df'),
                'function_heatmap_df': heatmap_data.get('func_heatmap_df'),
                'filtered_heatmap_df': heatmap_data.get('filtered_heatmap_df'),
                # Row/legend label. When more than one protein is being compared
                # the bare sample name ("Bitter") repeats across proteins and can't
                # tell them apart, so prefix the protein name ("Beta-casein Bitter");
                # a single protein keeps the sample name alone. `var_label` always
                # carries the bare sample name for callers (e.g. the differential
                # comparison track) that build their own protein-aware labels.
                'label': f"{protein_name} {var_key}" if len(selected_proteins) > 1 else var_key,
                'var_label': var_key,
                'is_func_df_all_none': (
                    heatmap_data.get('func_heatmap_df') is None
                    or (isinstance(heatmap_data.get('func_heatmap_df'), pd.DataFrame)
                        and heatmap_data['func_heatmap_df'].isnull().all().all())
                ),
            }

    return available, messages


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# The x-axis-label / short-protein-name logic is defined once in
# heatmap_renderer (`_derive_xaxis_label` / `_derive_protein_short_name`) and
# imported where needed — see generate_heatmap below.


# Generate heatmap plot
# ---------------------------------------------------------------------------

def generate_heatmap(
    available_data_variables_dict: dict,
    plot_params: dict,
) -> tuple:
    """
    Call heatmap_renderer.update_plot() and return
    (portrait_b64, landscape_b64, compact_json, messages).

    portrait_b64 / landscape_b64 – base64-encoded PNG strings (or None).
    compact_json                  – Plotly figure JSON string for the compact
                                    interactive plot (or None).
    messages                      – list of error / notification strings.
    """
    import base64

    try:
        from .heatmap_renderer import update_plot, _derive_xaxis_label
    except ImportError:
        return None, None, None, ['Could not import heatmap_renderer.']

    if not available_data_variables_dict:
        return None, None, None, ['No data available for plotting.']

    pp = dict(plot_params)

    # "Strip protein name" toggle (Appearance Settings), default on. When the
    # user unchecks it, swap the short display name for the full raw name so the
    # axis label, title, and comparison labels all show the unstripped form.
    if not pp.get('strip_protein_name', True):
        for vd in available_data_variables_dict.values():
            raw = vd.get('protein_name_raw')
            if raw:
                vd['protein_name'] = raw

    # Auto-derive the x-axis label from protein name(s) when the user has left
    # the field blank — mirrors the notebook's `protein_name_short + " Sequence"` logic.
    user_xaxis_label = pp.get('xaxis_label', '').strip()
    resolved_xaxis_label = user_xaxis_label if user_xaxis_label else _derive_xaxis_label(available_data_variables_dict)

    try:
        fig_port, fig_land, fig_compact, fig_port_plotly, fig_land_plotly, errors, notifications = update_plot(
            available_data_variables_dict,
            ms_average_choice=pp.get('ms_average_choice', 'yes'),
            bio_or_pep=pp.get('bio_or_pep', 'no'),
            selected_peptides=pp.get('selected_peptides', []),
            selected_functions=pp.get('selected_functions', []),
            hm_selected_color=pp.get('hm_selected_color', 'RdYlGn_r'),
            lp_selected_color=pp.get('lp_selected_color', 'Set3'),
            avglp_selected_color=pp.get('avglp_selected_color', 'Dark2'),
            xaxis_label=resolved_xaxis_label,
            # Pass the raw user value (blank = auto). The renderer picks the
            # fallback: effect-size metric name in comparison mode, else the
            # abundance default — a user label always overrides both.
            yaxis_label=(pp.get('yaxis_label') or '').strip(),
            yaxis_position=pp.get('yaxis_position', 5),
            legend_title_input_1=pp.get('legend_title_1', 'Sample Type:'),
            legend_title_input_2=pp.get('legend_title_2', 'Peptide Counts:'),
            legend_title_input_3=pp.get('legend_title_3', 'Abundance:'),
            # 'all-peptides' is the value the UI sends; 'All' was never a recognised
            # filter_type and fell through every branch in process_available_data.
            filter_type=pp.get('filter_type', 'all-peptides'),
            log_transform=pp.get('log_transform', False),
            manual_y_axis=pp.get('manual_y_axis', False),
            y_min_manual=pp.get('y_min', 0.0),
            y_max_manual=pp.get('y_max', 1.0),
            plot_compact=pp.get('plot_compact', False),
            plot_landscape_interactive=pp.get('plot_landscape_interactive', False),
            plot_portrait_interactive=pp.get('plot_portrait_interactive', False),
            # Optional user figure title; blank → renderers draw no title.
            plot_title=pp.get('plot_title', '').strip(),
            # Print the AA letter on each sample's colored tiles (portrait/landscape
            # Plotly) instead of the shared grey sequence strip.
            aa_on_tiles=pp.get('aa_on_tiles', False),
            # Differential-comparison track (R2-2d): when comparison_mode is on,
            # the row-1 panel shows a signed SMD/log2FC line for series_a vs
            # series_b instead of the abundance line.
            comparison_mode=pp.get('comparison_mode', False),
            series_a=pp.get('series_a'),
            series_b=pp.get('series_b'),
            comparison_metric=pp.get('comparison_metric', 'smd'),
            # Appearance Settings font-size overrides.
            font_size_xaxis_label=_coerce_font_size(pp.get('font_size_xaxis_label'), 14),
            font_size_yaxis_label=_coerce_font_size(pp.get('font_size_yaxis_label'), 14),
            font_size_legend=_coerce_font_size(pp.get('font_size_legend'), 14),
            font_size_plot_title=_coerce_font_size(pp.get('font_size_plot_title'), 16),
            font_size_xaxis_tick=_coerce_font_size(pp.get('font_size_xaxis_tick'), 14),
            font_size_yaxis_tick=_coerce_font_size(pp.get('font_size_yaxis_tick'), 14),
            font_size_var_label=_coerce_font_size(pp.get('font_size_var_label'), 14),
            # Legend placement: 'right' (one combined legend beside the plot) or
            # 'below' (one unit per legend group under the x-axis, side by side),
            # which frees the figure width for the heatmap itself.
            legend_position=pp.get('legend_position', 'right'),
        )
    except Exception as exc:
        return None, None, None, None, None, [f'Error generating heatmap: {exc}\n{traceback.format_exc()}']

    messages = list(errors or []) + list(notifications or [])

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import io as _io

    # Close all open figures before encoding to prevent memory accumulation.
    plt.close('all')

    def fig_to_b64(fig) -> str | None:
        if fig is None:
            return None
        buf = _io.BytesIO()
        # 300 dpi is the standard minimum for publication figures. The same
        # encoded image backs both the on-screen preview and the download, so
        # what a user exports matches what they see.
        fig.savefig(buf, format='png', dpi=EXPORT_DPI, bbox_inches='tight')
        buf.seek(0)
        data = base64.b64encode(buf.read()).decode()
        plt.close(fig)
        return data

    portrait_b64  = fig_to_b64(fig_port) if fig_port  is not None else None
    landscape_b64 = fig_to_b64(fig_land) if fig_land  is not None else None

    # Convert Plotly figures to JSON for transport to the frontend
    def _plotly_to_json(fig, name: str) -> str | None:
        if fig is None:
            return None
        try:
            import plotly.io as pio
            return pio.to_json(fig)
        except Exception as exc:
            messages.append(f'Warning: could not serialise {name} interactive plot: {exc}')
            return None

    compact_json           = _plotly_to_json(fig_compact,     'compact')
    portrait_interactive   = _plotly_to_json(fig_port_plotly, 'portrait')
    landscape_interactive  = _plotly_to_json(fig_land_plotly, 'landscape')

    return portrait_b64, landscape_b64, compact_json, portrait_interactive, landscape_interactive, messages
