"""
Django views for the data transformation wizard.
"""
from __future__ import annotations

import json
import os
import uuid

from django.conf import settings
from django.core.cache import cache
from django.http import JsonResponse, HttpResponse, FileResponse
from django.shortcuts import render
from django.views.decorators.http import require_POST, require_GET
from celery.result import AsyncResult

from .forms import PeptidomicUploadForm, GroupUploadForm, FastaUploadForm
from utils.toolbox_helpers import handle_uploaded_file, clear_temp_directory
from utils.lazy_import import LazyModule

# pandas/numpy/scipy (via these services, and via .tasks -> services.blast_search)
# are heavy; the URLconf import chain loads this module on the very first
# request handled after boot (probe or otherwise), so they're deferred until
# a view actually uses them.
pd = LazyModule('pandas')
data_loader = LazyModule('data_transformation.services.data_loader')
blast_search = LazyModule('data_transformation.services.blast_search')
group_processing = LazyModule('data_transformation.services.group_processing')
protein_handler = LazyModule('data_transformation.services.protein_handler')
data_combiner = LazyModule('data_transformation.services.data_combiner')
export_manager = LazyModule('data_transformation.services.export_manager')
tasks = LazyModule('data_transformation.tasks')


def _get_work_dir(request):
    """Get or create the work directory for this session."""
    work_dir = request.session.get('dt_work_dir')
    if work_dir and os.path.isdir(work_dir):
        return work_dir
    work_dir = blast_search.create_work_directory()
    # Make work dir accessible to celery_user (runs as different user than gunicorn)
    os.chmod(work_dir, 0o777)
    request.session['dt_work_dir'] = work_dir
    request.session['dt_session_id'] = str(uuid.uuid4())
    return work_dir


def _save_df(work_dir, name, df):
    """Save a DataFrame as pickle in the work directory."""
    path = os.path.join(work_dir, f'{name}.pkl')
    df.to_pickle(path)
    return path


def _load_df(work_dir, name):
    """Load a DataFrame from pickle in the work directory."""
    path = os.path.join(work_dir, f'{name}.pkl')
    if os.path.exists(path):
        return pd.read_pickle(path)
    return None


def _save_json(work_dir, name, data):
    """Save data as JSON in the work directory."""
    path = os.path.join(work_dir, f'{name}.json')
    with open(path, 'w') as f:
        json.dump(data, f)


def _load_json(work_dir, name):
    """Load data from JSON in the work directory."""
    path = os.path.join(work_dir, f'{name}.json')
    if os.path.exists(path):
        with open(path, 'r') as f:
            return json.load(f)
    return None


def _delete_json(work_dir, name):
    """Remove a saved JSON file from the work directory if it exists."""
    path = os.path.join(work_dir, f'{name}.json')
    if os.path.exists(path):
        os.remove(path)


def _parse_log_base(request):
    """Log base for correlation transforms: 10 (default) or 2. Log10 and Log2
    are mutually exclusive in the UI."""
    return 2 if str(request.GET.get('log_base', '10')).strip() == '2' else 10


# State derived from a specific peptidomic dataset. A fresh peptidomic upload
# invalidates all of it, so upload_files() clears these before repopulating —
# otherwise stale column renames / groups / protein mappings from a previous
# dataset are validated against the new one and fail ("columns not present").
_DOWNSTREAM_JSON = [
    'base_columns', 'raw_columns', 'column_renames', 'group_data',
    'tech_dup_mapping', 'annotation_source', 'threshold',
    'protein_decisions', 'source_renames', 'unresolvable_proteins',
]
_DOWNSTREAM_DF = [
    'functional_data', 'mbpdb_results', 'merged_df', 'pd_results_cleaned',
]


def _reset_downstream_state(request, work_dir):
    """Clear all session state derived from a previous peptidomic dataset so a
    new upload repopulates groups, renames, and protein mappings from scratch."""
    for name in _DOWNSTREAM_JSON:
        _delete_json(work_dir, name)
    for name in _DOWNSTREAM_DF:
        path = os.path.join(work_dir, f'{name}.pkl')
        if os.path.exists(path):
            os.remove(path)
    stale = os.path.join(work_dir, 'functional_file.orig')
    if os.path.exists(stale):
        os.remove(stale)
    for key in ('dt_blast_task_id', 'dt_uniprot_task_id'):
        request.session.pop(key, None)


def _apply_renames_to_list(cols, renames):
    """Apply a column-rename mapping to a list of column names.

    renames is stored as {new_name: original_column} (the user-facing JSON
    format). Any column matching an original is swapped for its new name;
    all others pass through unchanged.
    """
    if not renames:
        return list(cols)
    original_to_new = {orig: new for new, orig in renames.items()}
    return [original_to_new.get(c, c) for c in cols]


def _apply_renames_to_df(df, renames):
    """Rename a DataFrame's columns in place-safe fashion using {new: original}."""
    if df is None or not renames:
        return df
    original_to_new = {orig: new for new, orig in renames.items()}
    # Only rename columns that actually exist to avoid clobbering
    present = {orig: new for orig, new in original_to_new.items() if orig in df.columns}
    return df.rename(columns=present) if present else df


# ---------------------------------------------------------------------------
# Main wizard page
# ---------------------------------------------------------------------------

def wizard_view(request):
    """Render the single-page wizard template."""
    return render(request, 'peptide/data_transformation.html')


# ---------------------------------------------------------------------------
# Step 1: Upload & BLAST
# ---------------------------------------------------------------------------

@require_POST
def upload_files(request):
    """Handle file uploads and validate peptidomic data."""
    try:
        form = PeptidomicUploadForm(request.POST, request.FILES)
        if not form.is_valid():
            return JsonResponse({'error': '; '.join(
                e for errors in form.errors.values() for e in errors
            )}, status=400)

        work_dir = _get_work_dir(request)

        # Detect merge mode: multiple files under 'merge_files'
        merge_files = request.FILES.getlist('merge_files')
        merged_from = 0

        if merge_files and len(merge_files) >= 2:
            # --- Merge mode: load, validate, and merge multiple datasets ---
            loaded = []
            all_warnings = []
            for mf in merge_files:
                mf_df, mf_status, mf_error, mf_warn = data_loader.load_and_validate_file(
                    mf, mf.name, 'Peptidomic'
                )
                if mf_status == 'no':
                    return JsonResponse({
                        'error': f'File "{mf.name}": {mf_error}'
                    }, status=400)
                if mf_warn:
                    all_warnings.append(f'{mf.name}: {mf_warn}')
                loaded.append((mf.name, mf_df))

            df = data_loader.merge_peptidomic_datasets(loaded)
            merged_from = len(loaded)
            warning_msg = '; '.join(all_warnings) if all_warnings else ''
            pep_file_name = ' + '.join(mf.name for mf in merge_files)
        else:
            # --- Single file mode (original behaviour) ---
            pep_file = request.FILES.get('peptidomic_file')
            if not pep_file:
                return JsonResponse(
                    {'error': 'No peptidomic data file provided.'}, status=400
                )

            df, status, error_msg, warning_msg = data_loader.load_and_validate_file(
                pep_file, pep_file.name, 'Peptidomic'
            )
            if status == 'no':
                return JsonResponse({'error': error_msg}, status=400)
            pep_file_name = pep_file.name

        # A new peptidomic dataset invalidates everything derived from the
        # previous one — wipe it so groups / renames / protein mappings are
        # repopulated fresh below and in the later wizard steps.
        _reset_downstream_state(request, work_dir)

        _save_df(work_dir, 'pd_results', df)

        # Save column list for group processing
        columns = df.columns.tolist()
        _save_json(work_dir, 'columns', columns)

        # Functional annotation (file upload or database search) is handled in a
        # separate step — see start_blast_search. Here we only detect an embedded
        # 'function' column on re-uploaded transformed files (below).
        func_warnings = []
        has_mbpdb = False
        mbpdb_rows = 0

        # Load default protein dictionary from the headers-only reference file
        # (948 proteins: bovine, human, and other milk-relevant species).
        # protein_database.fasta (71 full-sequence proteins) is kept for BLAST/PEPEX.
        protein_dict = data_loader.parse_fasta_headers_file(settings.PROTEIN_HEADERS_FILE)

        # Re-apply a custom FASTA file previously uploaded in the Protein Mapping
        # step (see upload_fasta_file) so re-submitting this form doesn't lose it.
        prior_fasta_path = os.path.join(work_dir, 'fasta_file.orig')
        if os.path.exists(prior_fasta_path):
            with open(prior_fasta_path, 'rb') as fh:
                protein_dict.update(data_loader.parse_uploaded_fasta(fh.read()))

        _save_json(work_dir, 'protein_dict', protein_dict)

        # Extract sequences for BLAST
        sequences = data_loader.extract_sequences(df)

        # Filter out peptides with non-alphabetic characters before BLAST.
        # BLAST only accepts amino-acid letters; a single bad sequence makes
        # blastp fail and the whole batch silently returns zero matches.
        # We surface this as a warning so the user can still proceed (e.g. to
        # use the app without running the database search).
        invalid_peptides = data_loader.find_invalid_peptide_sequences(sequences)
        if invalid_peptides:
            invalid_set = set(invalid_peptides)
            sequences = [s for s in sequences if s not in invalid_set]
            preview = invalid_peptides[:10]
            preview_str = ', '.join(preview)
            extra = (
                f' (and {len(invalid_peptides) - len(preview)} more)'
                if len(invalid_peptides) > len(preview) else ''
            )
            func_warnings.append(
                f'{len(invalid_peptides)} peptide sequence(s) contain '
                f'non-alphabetic characters and will be excluded from the '
                f'database search: {preview_str}{extra}. The database search only '
                f'accepts sequences made of amino-acid letters (A–Z). To '
                f'include these peptides, remove the non-alphabetic '
                f'characters from your dataset and re-upload. You can still '
                f'continue without running the database search.'
            )

        # Auto-detect embedded function column in peptidomic data (re-uploaded transformed file)
        if not has_mbpdb and 'function' in df.columns and df['function'].notna().any():
            has_mbpdb = True
            mbpdb_rows = int(df['function'].notna().sum())

        all_warnings = [w for w in [warning_msg] + func_warnings if w]
        result = {
            'success': True,
            'rows': len(df),
            'columns': len(df.columns),
            'sequences': len(sequences),
            'has_mbpdb': has_mbpdb,
            'mbpdb_rows': mbpdb_rows,
            'warning': '; '.join(all_warnings) if all_warnings else None,
        }
        if merged_from:
            result['merged_from'] = merged_from

        # Store sequences for BLAST
        _save_json(work_dir, 'sequences', sequences)

        # Save original filenames so they can be restored on resume. The FASTA
        # filename is managed separately by upload_fasta_file (Protein Mapping
        # step) — preserve whatever is already on record for it here.
        merge_file_names = [mf.name for mf in merge_files] if merged_from else None
        existing_file_names = _load_json(work_dir, 'uploaded_file_names') or {}
        _save_json(work_dir, 'uploaded_file_names', {
            'peptidomic_file': pep_file_name if not merged_from else None,
            'merge_files': merge_file_names,
            'functional_file': None,  # cleared by _reset_downstream_state above
            'fasta_file': existing_file_names.get('fasta_file'),
        })

        # Save original file bytes so they can be reloaded on resume
        try:
            if not merge_files or len(merge_files) < 2:
                pep_file = request.FILES.get('peptidomic_file')
                if pep_file:
                    pep_file.seek(0)
                    with open(os.path.join(work_dir, 'peptidomic_file.orig'), 'wb') as fh:
                        fh.write(pep_file.read())
            else:
                # For merge mode, save the merged DataFrame as the "original"
                # so session resume works without re-uploading individual files
                _save_df(work_dir, 'pd_results', df)
                stale = os.path.join(work_dir, 'peptidomic_file.orig')
                if os.path.exists(stale):
                    os.remove(stale)

            # functional_file.orig is written by start_blast_search (Functional
            # Annotation step) — not touched here.
            # fasta_file.orig is managed separately by upload_fasta_file
            # (Protein Mapping step) — not touched here.
        except Exception:
            pass  # Non-fatal: resume will still show filenames without reloadable bytes

        return JsonResponse(result)

    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        # Return the actual exception details so the user can diagnose it
        return JsonResponse({
            'error': f'Upload failed: {type(e).__name__}: {str(e)}',
            'detail': tb.split('\n')[-3] if tb else '',
        }, status=500)


@require_POST
def start_blast_search(request):
    """Run the functional annotation step: either ingest an uploaded functional
    data file, or launch a database (BLAST) search as a Celery task."""
    work_dir = _get_work_dir(request)
    sequences = _load_json(work_dir, 'sequences')

    # Functional Annotation step inputs (sent with the "Search" action).
    func_file = request.FILES.get('functional_file')
    if func_file:
        func_df, f_status, f_error, _ = data_loader.load_and_validate_file(
            func_file, func_file.name, 'MBPDB'
        )
        if f_status == 'no':
            return JsonResponse(
                {'error': f'Functional data file: {f_error}'}, status=400
            )
        if func_df is not None:
            _save_df(work_dir, 'functional_data', func_df)
            file_names = _load_json(work_dir, 'uploaded_file_names') or {}
            file_names['functional_file'] = func_file.name
            _save_json(work_dir, 'uploaded_file_names', file_names)
            try:
                func_file.seek(0)
                with open(os.path.join(work_dir, 'functional_file.orig'), 'wb') as fh:
                    fh.write(func_file.read())
            except Exception:
                pass  # Non-fatal: resume still shows the filename

    posted_threshold = request.POST.get('similarity_threshold')
    if posted_threshold not in (None, ''):
        try:
            _save_json(work_dir, 'threshold', int(posted_threshold))
        except (TypeError, ValueError):
            pass

    source = request.POST.get('annotation_source')
    if source:
        _save_json(work_dir, 'annotation_source', source)

    threshold = _load_json(work_dir, 'threshold') or 80

    # Check for uploaded functional data (skip BLAST if provided)
    func_df = _load_df(work_dir, 'functional_data')
    if func_df is not None:
        _save_df(work_dir, 'mbpdb_results', func_df)
        return JsonResponse({
            'skipped': True,
            'message': 'Using uploaded functional data instead of a database search',
            'count': len(func_df),
        })

    # Check if peptidomic data has an embedded function column (re-uploaded transformed file)
    pd_df = _load_df(work_dir, 'pd_results')
    if pd_df is not None and 'function' in pd_df.columns and pd_df['function'].notna().any():
        id_col = next((c for c in ['Sequence', 'sequence'] if c in pd_df.columns), None)
        if id_col:
            synth = pd_df[[id_col, 'function']].dropna(subset=[id_col, 'function']).copy()
            synth = synth[synth['function'].astype(str).str.strip() != '']
            synth = synth.rename(columns={id_col: 'search_peptide'})
            synth = synth.drop_duplicates(subset=['search_peptide', 'function'])
            _save_df(work_dir, 'mbpdb_results', synth)
            count = len(synth)
        else:
            count = int(pd_df['function'].notna().sum())
        return JsonResponse({
            'skipped': True,
            'message': 'Using function data embedded in uploaded file',
            'count': count,
        })

    if not sequences:
        return JsonResponse({
            'error': 'No peptide sequences available. Upload peptidomic data first, '
                     'or provide a functional annotation file.'
        }, status=400)

    task = tasks.run_blast_search_task.delay(work_dir, sequences, threshold)
    request.session['dt_blast_task_id'] = task.id

    return JsonResponse({'task_id': task.id})


@require_GET
def get_blast_results(request, task_id):
    """Get BLAST search results after task completes."""
    task_result = AsyncResult(str(task_id))

    if not task_result.ready():
        return JsonResponse({'status': 'pending'})

    if task_result.failed():
        return JsonResponse({'status': 'failed', 'error': str(task_result.result)}, status=500)

    result_data = task_result.result
    work_dir = _get_work_dir(request)

    # Load result pickle
    result_path = result_data.get('result_path', '')
    if os.path.exists(result_path):
        mbpdb_df = pd.read_pickle(result_path)
        _save_df(work_dir, 'mbpdb_results', mbpdb_df)
        count = len(mbpdb_df)
    else:
        count = 0

    return JsonResponse({
        'status': 'complete',
        'count': count,
    })


# ---------------------------------------------------------------------------
# Step 2: Study Variable Grouping
# ---------------------------------------------------------------------------

@require_GET
def get_step2_form(request):
    """Return data needed for Step 2 (grouping)."""
    try:
        work_dir = _get_work_dir(request)
        columns = _load_json(work_dir, 'columns') or []

        # Filter columns for abundance selection
        df = _load_df(work_dir, 'pd_results')
        if df is not None:
            filtered = data_loader.get_filtered_columns(df)
            # Fallback: if filter is too aggressive, return all non-metadata columns
            if len(filtered) < 2:
                filtered = data_loader.get_filtered_columns_fallback(df)
        else:
            filtered = columns

        # base_columns: the original (pre-rename) abundance columns, saved once so
        # the rename picker always offers the true instrument names and renames can
        # be edited/removed without losing the source list.
        base_columns = _load_json(work_dir, 'base_columns')
        if not base_columns:
            base_columns = list(filtered)
            _save_json(work_dir, 'base_columns', base_columns)

        # Column renames: {new_name: original_column}. Layered on top of the base
        # so the tech-rep and grouping panels display the simplified names.
        column_renames = _load_json(work_dir, 'column_renames') or {}
        raw_columns = _apply_renames_to_list(base_columns, column_renames)
        _save_json(work_dir, 'raw_columns', raw_columns)

        # On resume: use the user's previously saved tech rep mapping instead of
        # re-running auto-detection (which would clobber any manual adjustments).
        # Auto-detection only runs on a clean first load with no renames applied.
        tech_dup_mapping = _load_json(work_dir, 'tech_dup_mapping')
        if tech_dup_mapping is None:
            tech_dup_mapping = {}
            if df is not None and filtered and not column_renames:
                _, tech_dup_mapping, filtered = group_processing.collapse_technical_duplicates(df, filtered)
                if tech_dup_mapping:
                    _save_json(work_dir, 'tech_dup_mapping', tech_dup_mapping)
            else:
                filtered = list(raw_columns)
        else:
            # Recompute post-collapse column list from saved mapping
            filtered = group_processing.compute_bio_rep_columns(raw_columns, tech_dup_mapping)

        # Return previously defined groups so the frontend can restore them.
        saved_group_data = _load_json(work_dir, 'group_data') or {}
        saved_groups = [
            {'name': info['grouping_variable'], 'columns': info['abundance_columns']}
            for info in saved_group_data.values()
            if info.get('grouping_variable') and info.get('abundance_columns')
        ]

        return JsonResponse({
            'columns': filtered,
            'all_columns': columns,
            'raw_columns': raw_columns,
            'base_columns': base_columns,
            'column_renames': column_renames,
            'tech_dup_mapping': tech_dup_mapping,
            'saved_groups': saved_groups,
        })
    except Exception as e:
        return JsonResponse({'error': f'Could not load step 2: {str(e)}', 'columns': [], 'all_columns': []}, status=500)


@require_POST
def submit_tech_reps(request):
    """
    Submit technical replicate groups (auto-detected and/or manually defined).

    The submitted groups are the source of truth — they fully replace any previously
    saved mapping. Auto-detected groups are pre-populated on the frontend so the user
    can review, remove, or add to them before submitting.
    """
    work_dir = _get_work_dir(request)

    try:
        body = json.loads(request.body)
    except json.JSONDecodeError:
        return JsonResponse({'error': 'Invalid JSON'}, status=400)

    submitted_groups = body.get('tech_reps', [])
    raw_columns = _load_json(work_dir, 'raw_columns') or []

    # Build mapping entirely from submitted groups (name non-empty, at least 2 columns)
    tech_dup_mapping = {}
    for group in submitted_groups:
        name = group.get('name', '').strip()
        cols = group.get('columns', [])
        if name and len(cols) >= 2:
            tech_dup_mapping[name] = cols

    _save_json(work_dir, 'tech_dup_mapping', tech_dup_mapping)

    # Compute the post-collapse column list for Phase 2
    bio_rep_columns = group_processing.compute_bio_rep_columns(raw_columns, tech_dup_mapping)

    return JsonResponse({
        'success': True,
        'columns': bio_rep_columns,
        'tech_dup_mapping': tech_dup_mapping,
    })


@require_POST
def upload_tech_rep_json(request):
    """Upload and parse a technical replicate definition JSON file."""
    if 'tech_rep_file' not in request.FILES:
        return JsonResponse({'error': 'No file provided'}, status=400)

    try:
        content = request.FILES['tech_rep_file'].read()
        data = json.loads(content)
    except (json.JSONDecodeError, UnicodeDecodeError):
        return JsonResponse({'error': 'Invalid JSON file'}, status=400)

    if not isinstance(data, dict):
        return JsonResponse({'error': 'File must be a JSON object mapping bio replicate names to lists of column names'}, status=400)

    groups = []
    errors = []
    for name, cols in data.items():
        if not isinstance(cols, list) or len(cols) < 2:
            errors.append(f'"{name}" must map to a list of at least 2 column names')
            continue
        groups.append({'name': str(name), 'columns': [str(c) for c in cols]})

    if errors:
        return JsonResponse({'error': '; '.join(errors)}, status=400)

    return JsonResponse({'success': True, 'groups': groups})


def _parse_rename_mapping(data, base_columns):
    """Validate a {new_name: original_column} mapping against base columns.

    Returns (renames, error). renames is a cleaned {new_name: original_column}
    dict; error is a message string or None.
    """
    if not isinstance(data, dict):
        return None, ('File must be a JSON object mapping new names to original '
                      'column names, e.g. {"T_0.2a": "Abundance F1 Sample Threshold"}')
    base_set = set(base_columns or [])
    renames = {}
    seen_originals = {}
    errors = []
    for new_name, original in data.items():
        new_name = str(new_name).strip()
        original = str(original).strip()
        if not new_name or not original:
            errors.append('New names and original columns must be non-empty')
            continue
        if base_set and original not in base_set:
            errors.append(f'"{original}" is not one of the available columns')
            continue
        if original in seen_originals:
            errors.append(f'"{original}" is mapped more than once')
            continue
        seen_originals[original] = new_name
        renames[new_name] = original
    if errors:
        return None, '; '.join(errors)
    return renames, None


@require_POST
def upload_rename_json(request):
    """Upload and parse a column-rename definition JSON file.

    Expected format: {new_name: original_column}, e.g.
    {"T_0.2a": "Abundance F1 Sample Threshold"}. Validated against the original
    (pre-rename) abundance columns but not persisted here — the frontend submits
    the final mapping via submit-renames/.
    """
    work_dir = _get_work_dir(request)
    if 'rename_file' not in request.FILES:
        return JsonResponse({'error': 'No file provided'}, status=400)

    try:
        data = json.loads(request.FILES['rename_file'].read())
    except (json.JSONDecodeError, UnicodeDecodeError):
        return JsonResponse({'error': 'Invalid JSON file'}, status=400)

    base_columns = _load_json(work_dir, 'base_columns') or []
    renames, error = _parse_rename_mapping(data, base_columns)
    if error:
        return JsonResponse({'error': error}, status=400)

    return JsonResponse({'success': True, 'renames': renames})


@require_POST
def submit_renames(request):
    """Persist the column-rename mapping ({new_name: original_column})."""
    work_dir = _get_work_dir(request)

    try:
        body = json.loads(request.body)
    except json.JSONDecodeError:
        return JsonResponse({'error': 'Invalid JSON'}, status=400)

    base_columns = _load_json(work_dir, 'base_columns') or []
    renames, error = _parse_rename_mapping(body.get('renames', {}), base_columns)
    if error:
        return JsonResponse({'error': error}, status=400)

    _save_json(work_dir, 'column_renames', renames)

    # Keep the effective raw column list in sync for the tech-rep step.
    raw_columns = _apply_renames_to_list(base_columns, renames)
    _save_json(work_dir, 'raw_columns', raw_columns)

    return JsonResponse({'success': True, 'renames': renames, 'raw_columns': raw_columns})


@require_GET
def download_rename_map(request):
    """Download the column-rename mapping as a reusable JSON key.

    Returns the saved {new_name: original_column} mapping, or an identity
    template built from the current columns if none has been defined yet.
    """
    work_dir = _get_work_dir(request)
    renames = _load_json(work_dir, 'column_renames')
    if not renames:
        base_columns = _load_json(work_dir, 'base_columns') or []
        renames = {col: col for col in base_columns}

    content = json.dumps(renames, indent=2)
    response = HttpResponse(content, content_type='application/json')
    response['Content-Disposition'] = 'attachment; filename="column_rename_key.json"'
    return response


@require_POST
def upload_group_json(request):
    """Upload and parse a group definition JSON file."""
    work_dir = _get_work_dir(request)
    columns = _load_json(work_dir, 'columns') or []
    # Validate uploaded groups against the renamed columns the user sees in the UI.
    column_renames = _load_json(work_dir, 'column_renames') or {}
    columns = _apply_renames_to_list(columns, column_renames)

    if 'group_file' not in request.FILES:
        return JsonResponse({'error': 'No file provided'}, status=400)

    content = request.FILES['group_file'].read()
    group_data, error = group_processing.parse_group_json(content, columns)

    if error:
        # A common mistake: loading the group-definition example (which uses
        # simplified sample names) before uploading/defining the column rename
        # key, so none of those names exist in the raw instrument columns yet.
        if not column_renames and 'not present in the current dataset' in error:
            error += (
                '\n\nIt appears you uploaded the group definitions with simplified / '
                'renamed columns before renaming the columns in step 1, "Rename '
                'Abundance Columns", above. Define your column renames there first, '
                'then load the group example again.'
            )
        return JsonResponse({'error': error}, status=400)

    _save_json(work_dir, 'group_data', group_data)

    # Return name→columns mapping so JS can populate definedGroups correctly
    return JsonResponse({
        'success': True,
        'groups': {
            info['grouping_variable']: info['abundance_columns']
            for info in group_data.values()
        },
    })


@require_POST
def submit_groups(request):
    """Submit manually defined groups."""
    work_dir = _get_work_dir(request)

    try:
        body = json.loads(request.body)
    except json.JSONDecodeError:
        return JsonResponse({'error': 'Invalid JSON'}, status=400)

    groups = body.get('groups', [])
    if not groups:
        return JsonResponse({'error': 'No groups provided'}, status=400)

    group_data = {}
    for group_def in groups:
        group_data, error = group_processing.build_group_data(
            group_def.get('columns', []),
            group_def.get('name', ''),
            group_data
        )
        if error:
            return JsonResponse({'error': error}, status=400)

    _save_json(work_dir, 'group_data', group_data)

    return JsonResponse({
        'success': True,
        'groups': {
            gid: info['grouping_variable']
            for gid, info in group_data.items()
        },
    })


@require_POST
def reset_groups(request):
    """Clear any saved group definitions (uploaded or manually built).

    The Reset Groups button clears the client-side state, but the grouping
    variable is also persisted server-side in group_data.json. Without deleting
    it, a page resume / step-2 re-fetch restores the prior groups.
    """
    work_dir = _get_work_dir(request)
    _delete_json(work_dir, 'group_data')
    return JsonResponse({'success': True})


@require_POST
def skip_groups(request):
    """Skip grouping - create individual groups for selected columns."""
    work_dir = _get_work_dir(request)

    try:
        body = json.loads(request.body)
    except json.JSONDecodeError:
        return JsonResponse({'error': 'Invalid JSON'}, status=400)

    selected = body.get('columns', [])
    if not selected:
        # Use all filtered columns
        df = _load_df(work_dir, 'pd_results')
        if df is not None:
            selected = data_loader.get_filtered_columns(df)

    group_data, error = group_processing.build_no_group_data(selected)
    if error:
        return JsonResponse({'error': error}, status=400)

    _save_json(work_dir, 'group_data', group_data)

    return JsonResponse({'success': True, 'message': 'Individual column groups created'})


# ---------------------------------------------------------------------------
# Step 3: Protein Mapping
# ---------------------------------------------------------------------------

@require_GET
def get_step3_form(request):
    """Return data needed for Step 3 (protein mapping)."""
    try:
        work_dir = _get_work_dir(request)
        df = _load_df(work_dir, 'pd_results')
        protein_dict = _load_json(work_dir, 'protein_dict') or {}

        if df is None:
            return JsonResponse({'error': 'No data loaded'}, status=400)

        all_missing_ids = list(data_loader.find_missing_proteins(df, protein_dict))

        # Split missing IDs into those still fetchable and those already tried
        # but not found (non-standard format or absent from UniProt).
        unresolvable = set(_load_json(work_dir, 'unresolvable_proteins') or [])
        missing_ids = [pid for pid in all_missing_ids if pid not in unresolvable]
        unresolvable_ids = [pid for pid in all_missing_ids if pid in unresolvable]

        # Summarize known proteins found in the data
        all_protein_ids = set()
        if 'Protein' in df.columns:
            for val in df['Protein'].dropna():
                pid = data_loader.extract_protein_id(val)
                if isinstance(pid, list):
                    all_protein_ids.update(pid)
                elif pid:
                    all_protein_ids.add(pid)
        known_proteins = [
            {
                'id': pid,
                'name': protein_dict[pid].get('name', ''),
                'species': protein_dict[pid].get('species', ''),
            }
            for pid in sorted(all_protein_ids)
            if pid in protein_dict
        ]

        # Get protein combinations
        combinations = protein_handler.get_protein_combinations(df, protein_dict)
        combo_details = protein_handler.get_combination_details(
            combinations, df, protein_dict
        ) if combinations else []

        return JsonResponse({
            'missing_protein_count': len(missing_ids),
            'missing_ids': missing_ids[:50],
            'unresolvable_count': len(unresolvable_ids),
            'unresolvable_ids': unresolvable_ids[:50],
            'known_protein_count': len(known_proteins),
            'known_proteins': known_proteins[:50],
            'combinations': combo_details,
            'has_combinations': len(combo_details) > 0,
            'saved_decisions': _decisions_to_ui(_load_json(work_dir, 'protein_decisions') or {}),
            'protein_sources': protein_handler.get_protein_sources(df, protein_dict),
            'saved_source_renames': _load_json(work_dir, 'source_renames') or [],
        })
    except Exception as e:
        return JsonResponse({'error': f'Could not load step 3: {str(e)}'}, status=500)


@require_POST
def upload_fasta_file(request):
    """
    Upload and validate a custom FASTA file, merging its protein entries into
    protein_dict on top of whatever is already known (default headers file
    and/or previously fetched UniProt results).
    """
    work_dir = _get_work_dir(request)
    form = FastaUploadForm(request.POST, request.FILES)
    if not form.is_valid():
        return JsonResponse({'error': '; '.join(
            e for errors in form.errors.values() for e in errors
        )}, status=400)

    fasta_file = request.FILES['fasta_file']
    content = fasta_file.read()

    is_valid, error_msg = data_loader.validate_fasta_format(content)
    if not is_valid:
        return JsonResponse({'error': f'Invalid FASTA file: {error_msg}'}, status=400)

    user_proteins = data_loader.parse_uploaded_fasta(content, fasta_file.name)
    if not user_proteins:
        return JsonResponse(
            {'error': 'No protein sequences could be parsed from this FASTA file.'}, status=400
        )

    protein_dict = _load_json(work_dir, 'protein_dict') or {}
    protein_dict.update(user_proteins)
    _save_json(work_dir, 'protein_dict', protein_dict)

    # Save the filename and original bytes so the upload survives session resume
    # and so re-submitting Step 1 can re-apply it (see upload_files).
    file_names = _load_json(work_dir, 'uploaded_file_names') or {}
    file_names['fasta_file'] = fasta_file.name
    _save_json(work_dir, 'uploaded_file_names', file_names)
    try:
        with open(os.path.join(work_dir, 'fasta_file.orig'), 'wb') as fh:
            fh.write(content)
    except Exception:
        pass  # Non-fatal: mapping is still applied, just won't survive resume

    return JsonResponse({'success': True, 'added': len(user_proteins)})


@require_POST
def start_uniprot_fetch(request):
    """Launch UniProt fetch as a Celery task."""
    work_dir = _get_work_dir(request)
    df = _load_df(work_dir, 'pd_results')
    protein_dict = _load_json(work_dir, 'protein_dict') or {}

    if df is None:
        return JsonResponse({'error': 'No data loaded'}, status=400)

    missing_ids = list(data_loader.find_missing_proteins(df, protein_dict))

    if not missing_ids:
        return JsonResponse({'skipped': True, 'message': 'No missing proteins to fetch'})

    task = tasks.fetch_uniprot_task.delay(missing_ids)
    request.session['dt_uniprot_task_id'] = task.id

    return JsonResponse({'task_id': task.id, 'count': len(missing_ids)})


@require_POST
def save_uniprot_results(request):
    """
    Persist completed UniProt task results into protein_dict on disk.
    Called by JS once the UniProt poll finishes, before re-rendering step 3.
    The client passes the task_id explicitly so we never rely on session state.
    """
    work_dir = _get_work_dir(request)

    # Prefer task_id from the request body; fall back to session for compatibility.
    try:
        body = json.loads(request.body)
        uniprot_task_id = body.get('task_id')
    except (json.JSONDecodeError, AttributeError):
        uniprot_task_id = None

    if not uniprot_task_id:
        uniprot_task_id = request.session.get('dt_uniprot_task_id')

    if not uniprot_task_id:
        return JsonResponse({'saved': 0})

    task_result = AsyncResult(str(uniprot_task_id))
    if task_result.ready() and not task_result.failed():
        raw = task_result.result or {}
    else:
        # Celery backend unavailable — fall back to the Django cache copy stored
        # by fetch_uniprot_task immediately before it returned.
        cached_found = cache.get(f'uniprot_found_{uniprot_task_id}', {})
        raw = {'found': cached_found} if cached_found else {}

    if not raw:
        return JsonResponse({'saved': 0})

    # Support both the new structured format {found, not_found, skipped}
    # and the legacy flat format {pid: info} for backward compatibility.
    if 'found' in raw:
        uniprot_results = raw['found']
        unresolvable_new = raw.get('not_found', []) + raw.get('skipped', [])
    else:
        uniprot_results = raw
        unresolvable_new = []

    protein_dict = _load_json(work_dir, 'protein_dict') or {}
    saved = 0
    saved_proteins = []
    for pid, info in uniprot_results.items():
        if pid not in protein_dict or not protein_dict[pid].get('name'):
            protein_dict[pid] = info
            saved += 1
            saved_proteins.append({
                'id': pid,
                'name': info.get('name', ''),
                'species': info.get('species', ''),
            })

    _save_json(work_dir, 'protein_dict', protein_dict)

    # Persist IDs that can never be resolved so the UI can distinguish them
    # from proteins that simply haven't been fetched yet.
    if unresolvable_new:
        existing = _load_json(work_dir, 'unresolvable_proteins') or []
        merged = list(set(existing + unresolvable_new))
        _save_json(work_dir, 'unresolvable_proteins', merged)

    return JsonResponse({'saved': saved, 'proteins': saved_proteins})


def _translate_decisions(raw_decisions):
    """
    Translate simplified JS-format decisions to backend format.
    {combo: {action, protein_ids?, protein_id?}} → {combo: {protein_id: decision_str}}
    Handles both the simplified format and legacy backend format transparently.
    """
    decisions = {}
    for combo, decision_data in raw_decisions.items():
        if not isinstance(decision_data, dict) or 'action' not in decision_data:
            # Already in backend format (legacy) — pass through
            decisions[combo] = decision_data
            continue
        action = decision_data.get('action', 'ASIS').upper()
        protein_id = decision_data.get('protein_id', '').strip()
        all_proteins = [p.strip() for p in combo.split(';') if p.strip()]
        if action == 'ASIS' or not all_proteins:
            decisions[combo] = {p: 'ASIS' for p in all_proteins}
        elif action == 'SPLIT':
            selected = set(decision_data.get('protein_ids', []))
            decisions[combo] = {
                p: ('NEW' if p in selected else 'REMOVE') for p in all_proteins
            }
        elif action == 'NEW' and protein_id:
            decisions[combo] = {
                p: ('NEW' if p == protein_id else 'REMOVE') for p in all_proteins
            }
        elif action == 'CUSTOM' and protein_id:
            decisions[combo] = {all_proteins[0]: f'CUSTOM:{protein_id}'}
            for p in all_proteins[1:]:
                decisions[combo][p] = 'REMOVE'
        else:
            decisions[combo] = {p: 'ASIS' for p in all_proteins}
    return decisions


def _decisions_to_ui(raw_decisions):
    """Normalize saved protein decisions into UI restore-hints.

    Accepts either the simplified UI format ({combo: {action, protein_ids?,
    protein_id?}}, as saved by the Step-3 form) or the backend format
    ({combo: {protein_id: 'ASIS'|'NEW'|'REMOVE'|'CUSTOM:<id>'}}, as saved by an
    uploaded mapping key). Returns a format-agnostic dict the frontend can use to
    re-select the mode and repopulate the custom ID / split selection:

        {combo: {'mode': 'asis'|'split'|'custom',
                 'protein_ids': [...],      # for split
                 'protein_id': '<custom>'}} # for custom
    """
    ui = {}
    for combo, data in (raw_decisions or {}).items():
        if not isinstance(data, dict):
            continue
        if 'action' in data:
            # Simplified UI format
            action = str(data.get('action', 'ASIS')).upper()
            if action == 'SPLIT':
                ui[combo] = {'mode': 'split',
                             'protein_ids': list(data.get('protein_ids', []))}
            elif action == 'NEW' and data.get('protein_id'):
                ui[combo] = {'mode': 'split',
                             'protein_ids': [data['protein_id'].strip()]}
            elif action == 'CUSTOM' and data.get('protein_id'):
                ui[combo] = {'mode': 'custom',
                             'protein_id': data['protein_id'].strip()}
                if data.get('protein_name'):
                    ui[combo]['protein_name'] = str(data['protein_name']).strip()
            else:
                ui[combo] = {'mode': 'asis'}
        else:
            # Backend format: {protein_id: decision_str}
            custom_id = None
            split_ids = []
            for pid, dec in data.items():
                dec_str = str(dec)
                if dec_str.startswith('CUSTOM:'):
                    custom_id = dec_str.split('CUSTOM:', 1)[1].strip()
                elif dec_str.upper() == 'NEW':
                    split_ids.append(pid)
            if custom_id:
                ui[combo] = {'mode': 'custom', 'protein_id': custom_id}
            elif split_ids:
                ui[combo] = {'mode': 'split', 'protein_ids': split_ids}
            else:
                ui[combo] = {'mode': 'asis'}
    return ui


def _apply_custom_names(raw_decisions, protein_dict):
    """
    Persist any custom protein name entered alongside a CUSTOM decision into
    protein_dict, so ``add_protein_info`` picks it up at process time (it only
    overwrites 'Unknown' cells, and add_protein_info initializes protein_name
    to 'Unknown' before looking it up). Only the simplified UI format
    ({combo: {action, protein_id, protein_name}}) carries names; legacy
    backend-format decisions have no room for one.
    """
    for data in (raw_decisions or {}).values():
        if not isinstance(data, dict) or 'action' not in data:
            continue
        if str(data.get('action', '')).upper() != 'CUSTOM':
            continue
        protein_id = str(data.get('protein_id', '')).strip()
        protein_name = str(data.get('protein_name', '')).strip()
        if protein_id and protein_name:
            protein_dict.setdefault(protein_id, {})['name'] = protein_name


def _merge_uniprot_into_dict(request, protein_dict):
    """Merge any completed UniProt task results into protein_dict in place."""
    uniprot_task_id = request.session.get('dt_uniprot_task_id')
    if not uniprot_task_id:
        return
    task_result = AsyncResult(str(uniprot_task_id))
    if task_result.ready() and not task_result.failed():
        raw_uni = task_result.result or {}
    else:
        cached_found = cache.get(f'uniprot_found_{uniprot_task_id}', {})
        raw_uni = {'found': cached_found} if cached_found else {}
    uniprot_results = raw_uni.get('found', raw_uni) if raw_uni else {}
    for pid, info in uniprot_results.items():
        if pid not in protein_dict:
            protein_dict[pid] = info


@require_POST
def submit_protein_decisions(request):
    """Submit protein combination decisions."""
    work_dir = _get_work_dir(request)

    try:
        body = json.loads(request.body)
    except json.JSONDecodeError:
        return JsonResponse({'error': 'Invalid JSON'}, status=400)

    raw_decisions = body.get('decisions', {})
    source_renames = body.get('source_renames', []) or []
    protein_dict = _load_json(work_dir, 'protein_dict') or {}

    # Persist raw decisions + merge/rename groups so they can be downloaded as a
    # reusable mapping key and restored on return to this step.
    _save_json(work_dir, 'protein_decisions', raw_decisions)
    _save_json(work_dir, 'source_renames', source_renames)

    decisions = _translate_decisions(raw_decisions)

    _merge_uniprot_into_dict(request, protein_dict)
    _apply_custom_names(raw_decisions, protein_dict)

    df = _load_df(work_dir, 'pd_results')
    if df is None:
        return JsonResponse({'error': 'No data loaded'}, status=400)

    processed_df, error = protein_handler.apply_protein_decisions(
        df, decisions, protein_dict
    )
    if error:
        return JsonResponse({'error': error}, status=400)

    # Merge/rename accessions run last (authoritative) on the cleaned frame.
    if source_renames:
        processed_df, error = protein_handler.apply_source_renames(
            processed_df, source_renames, protein_dict
        )
        if error:
            return JsonResponse({'error': error}, status=400)

    _save_df(work_dir, 'pd_results_cleaned', processed_df)
    _save_json(work_dir, 'protein_dict', protein_dict)

    warnings = protein_handler.detect_rename_conflicts(source_renames, raw_decisions)
    return JsonResponse({'success': True, 'warnings': warnings})


@require_POST
def reset_proteins(request):
    """Clear any saved protein mapping decisions and merge/rename groups.

    Mirrors reset_groups: the Reset Mapping button clears the client-side
    combinations/merge-group state, but decisions are also persisted
    server-side in protein_decisions.json / source_renames.json. Without
    deleting them, a page resume / step-3 re-fetch restores the prior mapping.

    Also forgets a previously uploaded custom FASTA file (its saved bytes and
    filename record) so it isn't silently re-applied on the next resume or
    Step 1 re-submit.
    """
    work_dir = _get_work_dir(request)
    _delete_json(work_dir, 'protein_decisions')
    _delete_json(work_dir, 'source_renames')

    fasta_path = os.path.join(work_dir, 'fasta_file.orig')
    if os.path.exists(fasta_path):
        os.remove(fasta_path)
    file_names = _load_json(work_dir, 'uploaded_file_names') or {}
    if file_names.get('fasta_file'):
        file_names['fasta_file'] = None
        _save_json(work_dir, 'uploaded_file_names', file_names)

    return JsonResponse({'success': True})


@require_GET
def download_protein_map(request):
    """
    Download the protein mapping decisions as a reusable JSON key.
    If decisions have been submitted, returns those; otherwise generates a
    template from current combinations using their default decisions.

    Merge/rename groups are folded into ``protein_decisions`` as
    ``action: "MERGE"`` entries (see protein_handler.merge_groups_as_decisions)
    so the file has one unified shape instead of a separate array.
    """
    work_dir = _get_work_dir(request)

    saved = _load_json(work_dir, 'protein_decisions')
    saved_renames = _load_json(work_dir, 'source_renames') or []
    if saved:
        protein_decisions = dict(saved)
    else:
        # Generate a template from current combination defaults
        df = _load_df(work_dir, 'pd_results')
        protein_dict = _load_json(work_dir, 'protein_dict') or {}
        if df is None:
            return JsonResponse({'error': 'No data loaded'}, status=400)

        combinations = protein_handler.get_protein_combinations(df, protein_dict)
        combo_details = (
            protein_handler.get_combination_details(combinations, df, protein_dict)
            if combinations else []
        )
        protein_decisions = {}
        for c in combo_details:
            new_ids = [p['id'] for p in c['proteins'] if p.get('default_decision') == 'new']
            all_ids = [p['id'] for p in c['proteins']]
            if not new_ids or new_ids == all_ids:
                protein_decisions[c['combo']] = {'action': 'ASIS'}
            else:
                protein_decisions[c['combo']] = {'action': 'SPLIT', 'protein_ids': new_ids}

    protein_decisions.update(protein_handler.merge_groups_as_decisions(saved_renames))
    output_data = {'version': 3, 'protein_decisions': protein_decisions}

    content = json.dumps(output_data, indent=2)
    response = HttpResponse(content, content_type='application/json')
    response['Content-Disposition'] = 'attachment; filename="protein_mapping_key.json"'
    return response


@require_POST
def upload_protein_map(request):
    """
    Upload a previously downloaded protein mapping key JSON and apply it,
    replacing the need to manually configure combinations.

    Accepts the current unified format (merge/rename groups embedded in
    protein_decisions as ``action: "MERGE"`` entries), the older
    {version, protein_decisions, source_renames} format, and a flat
    {combo: decision} file, so previously downloaded keys keep working.
    """
    work_dir = _get_work_dir(request)

    if 'map_file' not in request.FILES:
        return JsonResponse({'error': 'No file provided'}, status=400)

    try:
        content = json.loads(request.FILES['map_file'].read())
    except (json.JSONDecodeError, UnicodeDecodeError):
        return JsonResponse({'error': 'Invalid JSON file'}, status=400)

    source_renames = []
    if 'protein_decisions' in content and isinstance(content['protein_decisions'], dict):
        raw_decisions = content['protein_decisions']
        if isinstance(content.get('source_renames'), list):
            source_renames = content['source_renames']
    elif isinstance(content, dict) and content:
        raw_decisions = content
    else:
        return JsonResponse({'error': 'Unrecognised mapping file format'}, status=400)

    # Pull any action:"MERGE" entries (current format) out into merge groups;
    # combine with a legacy top-level source_renames array if present.
    raw_decisions, embedded_groups = protein_handler.extract_merge_groups(raw_decisions)
    source_renames = source_renames + embedded_groups

    protein_dict = _load_json(work_dir, 'protein_dict') or {}
    _merge_uniprot_into_dict(request, protein_dict)
    _apply_custom_names(raw_decisions, protein_dict)

    _save_json(work_dir, 'protein_decisions', raw_decisions)
    _save_json(work_dir, 'source_renames', source_renames)
    decisions = _translate_decisions(raw_decisions)

    df = _load_df(work_dir, 'pd_results')
    if df is None:
        return JsonResponse({'error': 'No data loaded'}, status=400)

    processed_df, error = protein_handler.apply_protein_decisions(
        df, decisions, protein_dict
    )
    if error:
        return JsonResponse({'error': error}, status=400)

    # Merge/rename accessions run last (authoritative), mirroring submit.
    if source_renames:
        processed_df, error = protein_handler.apply_source_renames(
            processed_df, source_renames, protein_dict
        )
        if error:
            return JsonResponse({'error': error}, status=400)

    _save_df(work_dir, 'pd_results_cleaned', processed_df)
    _save_json(work_dir, 'protein_dict', protein_dict)

    return JsonResponse({
        'success': True,
        'applied': len(raw_decisions) + len(source_renames),
        'warnings': protein_handler.detect_rename_conflicts(source_renames, raw_decisions),
    })


@require_POST
def skip_protein_mapping(request):
    """Skip protein mapping step."""
    work_dir = _get_work_dir(request)

    protein_dict = _load_json(work_dir, 'protein_dict') or {}
    _merge_uniprot_into_dict(request, protein_dict)

    # Copy pd_results as pd_results_cleaned; still apply any previously-saved
    # merge/rename groups so skipping the combinations doesn't silently drop them.
    df = _load_df(work_dir, 'pd_results')
    if df is not None:
        source_renames = _load_json(work_dir, 'source_renames') or []
        if source_renames:
            df, error = protein_handler.apply_source_renames(df, source_renames, protein_dict)
            if error:
                return JsonResponse({'error': error}, status=400)
        _save_df(work_dir, 'pd_results_cleaned', df)

    _save_json(work_dir, 'protein_dict', protein_dict)

    return JsonResponse({'success': True})


# ---------------------------------------------------------------------------
# Step 4: Process & Export
# ---------------------------------------------------------------------------

@require_POST
def process_data(request):
    """Process all data and generate the merged dataset."""
    try:
        work_dir = _get_work_dir(request)

        pd_results = _load_df(work_dir, 'pd_results')
        pd_results_cleaned = _load_df(work_dir, 'pd_results_cleaned')
        mbpdb_results = _load_df(work_dir, 'mbpdb_results')
        group_data = _load_json(work_dir, 'group_data')
        protein_dict = _load_json(work_dir, 'protein_dict') or {}
        tech_dup_mapping = _load_json(work_dir, 'tech_dup_mapping')
        column_renames = _load_json(work_dir, 'column_renames') or {}

        if pd_results is None:
            return JsonResponse({'error': 'No peptidomic data loaded'}, status=400)

        # Apply column renames first so every downstream reference (tech reps,
        # groups, exports) resolves against the simplified names the user defined.
        if column_renames:
            pd_results = _apply_renames_to_df(pd_results, column_renames)
            if pd_results_cleaned is not None:
                pd_results_cleaned = _apply_renames_to_df(pd_results_cleaned, column_renames)

        # Apply technical duplicate collapsing so biological replicate columns exist
        # in the dataframe before group averaging runs.
        if tech_dup_mapping:
            pd_results = group_processing.apply_tech_dup_mapping(pd_results, tech_dup_mapping)
            if pd_results_cleaned is not None:
                pd_results_cleaned = group_processing.apply_tech_dup_mapping(
                    pd_results_cleaned, tech_dup_mapping
                )

        final_df = data_combiner.process_data(
            pd_results, pd_results_cleaned, mbpdb_results, group_data, protein_dict
        )

        if final_df is None:
            diag = []
            if pd_results is not None:
                diag.append(f'Input rows: {len(pd_results)}, columns: {list(pd_results.columns[:10])}')
            if group_data:
                for gid, ginfo in group_data.items():
                    missing = [c for c in ginfo.get('abundance_columns', [])
                               if c not in (pd_results.columns if pd_results is not None else [])]
                    if missing:
                        diag.append(f'Group "{ginfo["grouping_variable"]}" references missing columns: {missing}')
            return JsonResponse({
                'error': 'Data processing failed',
                'detail': '\n'.join(diag) if diag else 'process_data() returned None with no exception.',
            }, status=500)

        _save_df(work_dir, 'merged_df', final_df)

        # Determine available exports (bool() converts numpy.bool_ to JSON-safe Python bool)
        has_mbpdb = bool(mbpdb_results is not None and not mbpdb_results.empty)
        has_groups = bool(group_data is not None and len(group_data) > 0)
        has_function = bool('function' in final_df.columns and final_df['function'].notna().any())
        has_tech_rep_correlation = bool(
            tech_dup_mapping and
            any(len(v) >= 2 for v in tech_dup_mapping.values())
        )

        # Human-readable reason each unavailable export is greyed out, surfaced
        # as hover text on the disabled row in Step 4.
        export_reasons = {}
        no_groups_msg = (
            'No categorical groups were assigned in Step 2 (Study Variables), so '
            'per-group results cannot be computed. Go back to Step 2 and add at '
            'least one group.'
        )
        if not has_groups:
            for _k in ('group_definitions', 'sequence_list', 'summed_peptide',
                       'protein_analysis', 'group_correlation', 'replicate_correlation'):
                export_reasons[_k] = no_groups_msg
        if not has_mbpdb:
            export_reasons['mbpdb_results'] = (
                'No functional annotation data was provided or found. In Step 1, '
                'upload a functional data file or run the functional database search.'
            )
        if not (has_function and has_groups):
            export_reasons['summed_function'] = (
                no_groups_msg if (has_function and not has_groups) else
                'No functional annotation data is present in the processed dataset, '
                'so abundance cannot be summed by functional category. Provide '
                'functional data in Step 1.'
            )
        if not has_tech_rep_correlation:
            export_reasons['tech_rep_correlation'] = (
                'No technical replicate groups were defined in Step 2, so '
                'within-sample (technical replicate) correlations cannot be '
                'computed. Define technical replicates in the Step 2 '
                '"Technical Replicate Assignment" section.'
            )
        if not tech_dup_mapping:
            export_reasons['tech_rep_key'] = (
                'No technical replicate groups were defined in Step 2, so there '
                'is no technical replicate key to export.'
            )
        if not column_renames:
            export_reasons['column_rename_key'] = (
                'No abundance columns were renamed in Step 2, so there is no '
                'column rename key to export.'
            )

        return JsonResponse({
            'success': True,
            'rows': int(len(final_df)),
            'columns': int(len(final_df.columns)),
            'export_reasons': export_reasons,
            'exports': {
                'mbpdb_results': has_mbpdb,
                'group_definitions': has_groups,
                'merged_dataset': True,
                'sequence_list': has_groups,
                'summed_peptide': has_groups,
                'protein_analysis': has_groups,
                'summed_function': has_function and has_groups,
                'group_correlation': has_groups,
                'replicate_correlation': has_groups,
                'tech_rep_correlation': has_tech_rep_correlation,
                'tech_rep_key': bool(tech_dup_mapping),
                'column_rename_key': bool(column_renames),
                'protein_map': True,
            }
        })
    except Exception as exc:
        import traceback
        return JsonResponse(
            {'error': f'Unexpected error during processing: {exc}',
             'detail': traceback.format_exc()},
            status=500
        )


@require_GET
def view_export(request, export_type):
    """Return tabular data for in-app viewing as JSON."""
    import math
    import io as _io

    work_dir = _get_work_dir(request)
    merged_df = _load_df(work_dir, 'merged_df')
    group_data = _load_json(work_dir, 'group_data')
    protein_dict = _load_json(work_dir, 'protein_dict') or {}
    mbpdb_results = _load_df(work_dir, 'mbpdb_results')

    MAX_ROWS = 500

    def safe_val(v):
        import pandas as _pd
        if v is None:
            return None
        try:
            if _pd.isna(v):
                return None
        except (TypeError, ValueError):
            pass
        if hasattr(v, 'item'):
            v = v.item()
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            return None
        if isinstance(v, (list, tuple)):
            return str(v)
        return v

    def df_to_sheet(df, name):
        truncated = len(df) > MAX_ROWS
        rows = [[safe_val(v) for v in row] for row in df.head(MAX_ROWS).itertuples(index=False)]
        return {
            'name': name,
            'columns': [str(c) for c in df.columns],
            'rows': rows,
            'total_rows': len(df),
            'truncated': truncated,
        }

    try:
        sheets = None

        if export_type == 'mbpdb_results':
            if mbpdb_results is None or mbpdb_results.empty:
                return JsonResponse({'error': 'No functional annotation results available'}, status=404)
            sheets = [df_to_sheet(mbpdb_results, 'Functional Annotation')]

        elif export_type == 'group_definitions':
            if not group_data:
                return JsonResponse({'error': 'No group definitions'}, status=404)
            rows = [[info['grouping_variable'], ', '.join(info['abundance_columns'])]
                    for info in group_data.values()]
            sheets = [{'name': 'Group Definitions',
                       'columns': ['Group Name', 'Abundance Columns'],
                       'rows': rows, 'total_rows': len(rows), 'truncated': False}]

        elif export_type == 'merged_dataset':
            if merged_df is None:
                return JsonResponse({'error': 'No merged dataset'}, status=404)
            sheets = [df_to_sheet(merged_df, 'Merged Dataset')]

        elif export_type == 'sequence_list':
            if merged_df is None or not group_data:
                return JsonResponse({'error': 'No data available'}, status=404)
            result_df = export_manager.extract_sequences_by_name(merged_df, group_data)
            if result_df.empty:
                return JsonResponse({'error': 'No sequence data'}, status=404)
            sheets = [df_to_sheet(result_df, 'Sequence List')]

        elif export_type == 'summed_peptide':
            if merged_df is None or not group_data:
                return JsonResponse({'error': 'No data available'}, status=404)
            data = export_manager.summed_peptide_results(merged_df, group_data)
            if not data:
                return JsonResponse({'error': 'No summed peptide data'}, status=404)
            summary_rows, rep_rows = [], []
            for g, v in data.items():
                summary_rows.append([g, safe_val(v['total_Abundance']),
                                      safe_val(v['abundance_sem']),
                                      safe_val(v['unique_peptides']), safe_val(v['count_sem'])])
                ri = v['replicate_data']
                for i, rep in enumerate(ri['abundance_columns']):
                    ab = ri['replicate_abundances'][i] if i < len(ri['replicate_abundances']) else 0
                    ct = ri['replicate_counts'][i] if i < len(ri['replicate_counts']) else 0
                    rep_rows.append([g, rep, safe_val(ab), safe_val(ct)])
            sheets = [
                {'name': 'Summary',
                 'columns': ['Group', 'Total Abundance', 'Abundance SEM', 'Unique Peptides', 'Count SEM'],
                 'rows': summary_rows, 'total_rows': len(summary_rows), 'truncated': False},
                {'name': 'Replicate Details',
                 'columns': ['Group', 'Replicate', 'Total Abundance', 'Unique Peptides'],
                 'rows': rep_rows, 'total_rows': len(rep_rows), 'truncated': False},
            ]

        elif export_type == 'protein_analysis':
            if merged_df is None or not group_data:
                return JsonResponse({'error': 'No data available'}, status=404)
            xlsx_bytes, _ = export_manager.export_protein_data(merged_df, group_data, protein_dict)
            if xlsx_bytes is None:
                return JsonResponse({'error': 'No protein data'}, status=404)
            import openpyxl
            wb = openpyxl.load_workbook(_io.BytesIO(xlsx_bytes))
            sheets = []
            for sheet_name in wb.sheetnames:
                ws = wb[sheet_name]
                all_rows = list(ws.values)
                if not all_rows:
                    sheets.append({'name': sheet_name, 'columns': [], 'rows': [],
                                   'total_rows': 0, 'truncated': False})
                    continue
                cols = [str(c) if c is not None else '' for c in all_rows[0]]
                data_rows = [[safe_val(c) for c in row] for row in all_rows[1:]]
                sheets.append({'name': sheet_name, 'columns': cols,
                               'rows': data_rows[:MAX_ROWS],
                               'total_rows': len(data_rows),
                               'truncated': len(data_rows) > MAX_ROWS})

        elif export_type in ('summed_function', 'group_correlation', 'replicate_correlation'):
            if merged_df is None or not group_data:
                return JsonResponse({'error': 'No data available'}, status=404)
            correlation_type = request.GET.get('correlation_type', 'Pearson')
            log_transform = request.GET.get('log_transform', 'true').lower() == 'true'
            log_base = _parse_log_base(request)
            if export_type == 'summed_function':
                content = export_manager.export_summed_function_data(merged_df, group_data)
            elif export_type == 'group_correlation':
                content = export_manager.export_group_correlation(
                    merged_df, group_data, correlation_type, log_transform, log_base)
            else:
                content = export_manager.export_replicate_correlation(
                    merged_df, group_data, correlation_type, log_transform, log_base)
            if content is None:
                return JsonResponse({'error': 'No data available'}, status=404)
            import openpyxl
            wb = openpyxl.load_workbook(_io.BytesIO(content))
            sheets = []
            for sheet_name in wb.sheetnames:
                ws = wb[sheet_name]
                all_rows = list(ws.values)
                if not all_rows:
                    sheets.append({'name': sheet_name, 'columns': [], 'rows': [],
                                   'total_rows': 0, 'truncated': False})
                    continue
                cols = [str(c) if c is not None else '' for c in all_rows[0]]
                data_rows = [[safe_val(c) for c in row] for row in all_rows[1:]]
                sheets.append({'name': sheet_name, 'columns': cols,
                               'rows': data_rows[:MAX_ROWS],
                               'total_rows': len(data_rows),
                               'truncated': len(data_rows) > MAX_ROWS})

        elif export_type == 'tech_rep_correlation':
            pd_results = _load_df(work_dir, 'pd_results')
            tech_dup_mapping = _load_json(work_dir, 'tech_dup_mapping')
            if pd_results is None or not tech_dup_mapping:
                return JsonResponse({'error': 'No technical replicate data available'}, status=404)
            correlation_type = request.GET.get('correlation_type', 'Pearson')
            log_transform = request.GET.get('log_transform', 'true').lower() == 'true'
            content = export_manager.export_tech_rep_correlation(
                pd_results, tech_dup_mapping, correlation_type, log_transform,
                _parse_log_base(request))
            if content is None:
                return JsonResponse({'error': 'No valid technical replicate pairs'}, status=404)
            import openpyxl
            wb = openpyxl.load_workbook(_io.BytesIO(content))
            sheets = []
            for sheet_name in wb.sheetnames:
                ws = wb[sheet_name]
                all_rows = list(ws.values)
                if not all_rows:
                    sheets.append({'name': sheet_name, 'columns': [], 'rows': [],
                                   'total_rows': 0, 'truncated': False})
                    continue
                cols = [str(c) if c is not None else '' for c in all_rows[0]]
                data_rows = [[safe_val(c) for c in row] for row in all_rows[1:]]
                sheets.append({'name': sheet_name, 'columns': cols,
                               'rows': data_rows[:MAX_ROWS],
                               'total_rows': len(data_rows),
                               'truncated': len(data_rows) > MAX_ROWS})

        elif export_type == 'tech_rep_key':
            tech_dup_mapping = _load_json(work_dir, 'tech_dup_mapping')
            if not tech_dup_mapping:
                return JsonResponse({'error': 'No technical replicate mapping available'}, status=404)
            rows = [
                [bio_rep, ', '.join(tech_cols)]
                for bio_rep, tech_cols in tech_dup_mapping.items()
            ]
            sheets = [{'name': 'Technical Replicate Key',
                       'columns': ['Bio Replicate Name', 'Technical Replicate Columns'],
                       'rows': rows, 'total_rows': len(rows), 'truncated': False}]

        elif export_type == 'column_rename_key':
            column_renames = _load_json(work_dir, 'column_renames')
            if not column_renames:
                return JsonResponse({'error': 'No column renames defined'}, status=404)
            rows = [
                [new_name, original]
                for new_name, original in column_renames.items()
            ]
            sheets = [{'name': 'Column Rename Key',
                       'columns': ['New Name', 'Original Column'],
                       'rows': rows, 'total_rows': len(rows), 'truncated': False}]

        elif export_type == 'protein_map':
            saved = _load_json(work_dir, 'protein_decisions')
            saved_renames = _load_json(work_dir, 'source_renames') or []
            if saved:
                raw_mapping = dict(saved)
            else:
                # Generate template from current combinations (mirrors download_protein_map)
                df = _load_df(work_dir, 'pd_results')
                if df is None:
                    return JsonResponse({'error': 'No data loaded'}, status=404)
                combinations = protein_handler.get_protein_combinations(df, protein_dict)
                combo_details = (
                    protein_handler.get_combination_details(combinations, df, protein_dict)
                    if combinations else []
                )
                raw_mapping = {}
                for c in combo_details:
                    new_ids = [p['id'] for p in c['proteins'] if p.get('default_decision') == 'new']
                    all_ids = [p['id'] for p in c['proteins']]
                    if not new_ids or new_ids == all_ids:
                        raw_mapping[c['combo']] = {'action': 'ASIS'}
                    else:
                        raw_mapping[c['combo']] = {'action': 'SPLIT', 'protein_ids': new_ids}

            raw_mapping.update(protein_handler.merge_groups_as_decisions(saved_renames))

            rows = []
            for combo, decision in raw_mapping.items():
                if isinstance(decision, dict):
                    action = decision.get('action', 'ASIS')
                    if action == 'SPLIT':
                        details = ', '.join(decision.get('protein_ids', []))
                    elif action in ('NEW', 'CUSTOM'):
                        details = decision.get('protein_id', '')
                        if action == 'CUSTOM' and decision.get('protein_name'):
                            details += f" ({decision['protein_name']})"
                    elif action == 'MERGE':
                        sources = ', '.join(decision.get('protein_ids', []))
                        target = decision.get('target_id', '')
                        if decision.get('target_name'):
                            target += f" ({decision['target_name']})"
                        details = f"{sources} -> {target}"
                    else:
                        details = ''
                else:
                    action = str(decision)
                    details = ''
                rows.append([combo, action, details])

            sheets = [{'name': 'Protein Mapping Key',
                       'columns': ['Protein Combination', 'Action', 'Selected Protein(s)'],
                       'rows': rows, 'total_rows': len(rows), 'truncated': False}]

        if sheets is None:
            return JsonResponse({'error': 'Unknown export type'}, status=404)

        return JsonResponse({'sheets': sheets})

    except Exception as e:
        import traceback
        return JsonResponse({'error': f'View failed: {str(e)}',
                             'detail': traceback.format_exc().split('\n')[-3]}, status=500)


@require_GET
def download_export(request, export_type):
    """Download an export file."""
    work_dir = _get_work_dir(request)
    merged_df = _load_df(work_dir, 'merged_df')
    group_data = _load_json(work_dir, 'group_data')
    protein_dict = _load_json(work_dir, 'protein_dict') or {}
    mbpdb_results = _load_df(work_dir, 'mbpdb_results')

    correlation_type = request.GET.get('correlation_type', 'Pearson')
    log_transform = request.GET.get('log_transform', 'true').lower() == 'true'
    log_base = _parse_log_base(request)
    log_suffix = f'_log{log_base}' if log_transform else ''

    content = None
    filename = ''
    content_type = 'application/octet-stream'

    if export_type == 'mbpdb_results':
        content = export_manager.export_mbpdb_results(mbpdb_results)
        filename = 'functional_annotation_results.tsv'
        content_type = 'text/tab-separated-values'

    elif export_type == 'group_definitions':
        content = export_manager.export_group_definitions(group_data)
        filename = 'categorical_variable_definitions.json'
        content_type = 'application/json'

    elif export_type == 'merged_dataset':
        content = export_manager.export_merged_dataset(merged_df)
        filename = 'merged_dataframe.csv'
        content_type = 'text/csv'

    elif export_type == 'sequence_list':
        content = export_manager.export_sequence_list(merged_df, group_data)
        filename = 'list_of_peptides_by_sequences.csv'
        content_type = 'text/csv'

    elif export_type == 'summed_peptide':
        content = export_manager.export_summed_peptide_data(merged_df, group_data)
        filename = 'summed_peptide_results.xlsx'
        content_type = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'

    elif export_type == 'protein_analysis':
        content, _ = export_manager.export_protein_data(merged_df, group_data, protein_dict)
        filename = 'protein_analysis.xlsx'
        content_type = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'

    elif export_type == 'summed_function':
        content = export_manager.export_summed_function_data(merged_df, group_data)
        filename = 'functional_analysis.xlsx'
        content_type = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'

    elif export_type == 'group_correlation':
        content = export_manager.export_group_correlation(
            merged_df, group_data, correlation_type, log_transform, log_base
        )
        filename = f'group_correlations_{correlation_type.lower()}{log_suffix}.xlsx'
        content_type = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'

    elif export_type == 'replicate_correlation':
        content = export_manager.export_replicate_correlation(
            merged_df, group_data, correlation_type, log_transform, log_base
        )
        filename = f'replicate_correlations_{correlation_type.lower()}{log_suffix}.xlsx'
        content_type = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'

    elif export_type == 'tech_rep_key':
        tech_dup_mapping = _load_json(work_dir, 'tech_dup_mapping')
        if tech_dup_mapping:
            content = json.dumps(tech_dup_mapping, indent=2).encode('utf-8')
        filename = 'technical_replicate_key.json'
        content_type = 'application/json'

    elif export_type == 'column_rename_key':
        column_renames = _load_json(work_dir, 'column_renames')
        if column_renames:
            content = json.dumps(column_renames, indent=2).encode('utf-8')
        filename = 'column_rename_key.json'
        content_type = 'application/json'

    elif export_type == 'tech_rep_correlation':
        pd_results = _load_df(work_dir, 'pd_results')
        tech_dup_mapping = _load_json(work_dir, 'tech_dup_mapping')
        content = export_manager.export_tech_rep_correlation(
            pd_results, tech_dup_mapping, correlation_type, log_transform, log_base
        )
        filename = f'tech_rep_correlations_{correlation_type.lower()}{log_suffix}.xlsx'
        content_type = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'

    if content is None:
        return JsonResponse({'error': 'Export not available'}, status=404)

    response = HttpResponse(content, content_type=content_type)
    response['Content-Disposition'] = f'attachment; filename="{filename}"'
    return response


# ---------------------------------------------------------------------------
# Session state (for resume banner)
# ---------------------------------------------------------------------------

@require_GET
def session_state(request):
    """Return which steps have saved state so the wizard can offer to resume."""
    work_dir = request.session.get('dt_work_dir')
    if not work_dir or not os.path.isdir(work_dir):
        return JsonResponse({'has_session': False})
    file_names = _load_json(work_dir, 'uploaded_file_names') or {}
    return JsonResponse({
        'has_session': True,
        'has_data':      os.path.exists(os.path.join(work_dir, 'pd_results.pkl')),
        'has_groups':    os.path.exists(os.path.join(work_dir, 'group_data.json')),
        'has_processed': os.path.exists(os.path.join(work_dir, 'merged_df.pkl')),
        'file_names': file_names,
        'has_saved_files': {
            'peptidomic_file': os.path.exists(os.path.join(work_dir, 'peptidomic_file.orig')),
            'functional_file': os.path.exists(os.path.join(work_dir, 'functional_file.orig')),
            'fasta_file': os.path.exists(os.path.join(work_dir, 'fasta_file.orig')),
        },
    })


# ---------------------------------------------------------------------------
# Download All (ZIP)
# ---------------------------------------------------------------------------

@require_GET
def download_all_exports(request):
    """Download all available exports bundled as a single ZIP file."""
    import io
    import zipfile

    work_dir = _get_work_dir(request)
    merged_df = _load_df(work_dir, 'merged_df')
    if merged_df is None:
        return JsonResponse(
            {'error': 'No processed data available. Run Process Data first.'}, status=400
        )

    group_data = _load_json(work_dir, 'group_data')
    protein_dict = _load_json(work_dir, 'protein_dict') or {}
    mbpdb_results = _load_df(work_dir, 'mbpdb_results')
    tech_dup_mapping = _load_json(work_dir, 'tech_dup_mapping')
    pd_results = _load_df(work_dir, 'pd_results')

    # Each correlation export carries its own inline controls on the page, passed
    # here namespaced by export key. Fall back to the shared params for callers
    # (tests, scripts) that still send a single set.
    def _corr_opts(key):
        ct = (request.GET.get(f'{key}_correlation_type')
              or request.GET.get('correlation_type', 'Pearson'))
        lt_raw = request.GET.get(f'{key}_log_transform',
                                 request.GET.get('log_transform', 'true'))
        lt = str(lt_raw).lower() == 'true'
        lb_raw = request.GET.get(f'{key}_log_base', request.GET.get('log_base', '10'))
        lb = 2 if str(lb_raw).strip() == '2' else 10
        suffix = f'_log{lb}' if lt else ''
        return ct, lt, lb, suffix

    has_mbpdb = mbpdb_results is not None and not mbpdb_results.empty
    has_groups = bool(group_data and len(group_data) > 0)
    has_function = 'function' in merged_df.columns and merged_df['function'].notna().any()
    has_tech_rep = bool(tech_dup_mapping and any(len(v) >= 2 for v in tech_dup_mapping.values()))

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zf:
        def _add(filename, content):
            if content is not None:
                zf.writestr(filename, content)

        if has_mbpdb:
            _add('functional_annotation_results.tsv', export_manager.export_mbpdb_results(mbpdb_results))

        if has_groups:
            _add('categorical_variable_definitions.json',
                 export_manager.export_group_definitions(group_data))

        _add('merged_dataframe.csv', export_manager.export_merged_dataset(merged_df))

        if has_groups:
            _add('list_of_peptides_by_sequences.csv',
                 export_manager.export_sequence_list(merged_df, group_data))
            _add('summed_peptide_results.xlsx',
                 export_manager.export_summed_peptide_data(merged_df, group_data))
            protein_content, _ = export_manager.export_protein_data(merged_df, group_data, protein_dict)
            _add('protein_analysis.xlsx', protein_content)

            if has_function:
                _add('functional_analysis.xlsx',
                     export_manager.export_summed_function_data(merged_df, group_data))

            g_ct, g_lt, g_lb, g_sfx = _corr_opts('group_correlation')
            r_ct, r_lt, r_lb, r_sfx = _corr_opts('replicate_correlation')
            _add(f'group_correlations_{g_ct.lower()}{g_sfx}.xlsx',
                 export_manager.export_group_correlation(
                     merged_df, group_data, g_ct, g_lt, g_lb))
            _add(f'replicate_correlations_{r_ct.lower()}{r_sfx}.xlsx',
                 export_manager.export_replicate_correlation(
                     merged_df, group_data, r_ct, r_lt, r_lb))

        if has_tech_rep and pd_results is not None:
            t_ct, t_lt, t_lb, t_sfx = _corr_opts('tech_rep_correlation')
            _add(f'tech_rep_correlations_{t_ct.lower()}{t_sfx}.xlsx',
                 export_manager.export_tech_rep_correlation(
                     pd_results, tech_dup_mapping, t_ct, t_lt, t_lb))

        if tech_dup_mapping:
            _add('technical_replicate_key.json',
                 json.dumps(tech_dup_mapping, indent=2).encode('utf-8'))

        column_renames = _load_json(work_dir, 'column_renames')
        if column_renames:
            _add('column_rename_key.json',
                 json.dumps(column_renames, indent=2).encode('utf-8'))

        saved_decisions = dict(_load_json(work_dir, 'protein_decisions') or {})
        saved_renames = _load_json(work_dir, 'source_renames') or []
        saved_decisions.update(protein_handler.merge_groups_as_decisions(saved_renames))
        pm = {'version': 3, 'protein_decisions': saved_decisions}
        _add('protein_mapping_key.json', json.dumps(pm, indent=2).encode('utf-8'))

    buf.seek(0)
    response = HttpResponse(buf.getvalue(), content_type='application/zip')
    response['Content-Disposition'] = 'attachment; filename="data_transformation_exports.zip"'
    return response


# ---------------------------------------------------------------------------
# Reload saved files (serve previously uploaded raw bytes back to browser)
# ---------------------------------------------------------------------------

@require_GET
def reload_saved_file(request, file_key):
    """
    Return the raw bytes of a previously uploaded file so the browser can
    re-inject it into the file input and resubmit the upload form.
    file_key: 'peptidomic_file' | 'functional_file' | 'fasta_file'
    """
    allowed_keys = {'peptidomic_file', 'functional_file', 'fasta_file'}
    if file_key not in allowed_keys:
        return HttpResponse('Invalid file key.', status=400)

    work_dir = request.session.get('dt_work_dir')
    if not work_dir or not os.path.isdir(work_dir):
        return HttpResponse('No session found.', status=404)

    orig_path = os.path.join(work_dir, f'{file_key}.orig')
    if not os.path.exists(orig_path):
        return HttpResponse('Saved file not found.', status=404)

    file_names = _load_json(work_dir, 'uploaded_file_names') or {}
    original_name = file_names.get(file_key, f'{file_key}.bin')

    with open(orig_path, 'rb') as fh:
        data = fh.read()

    response = HttpResponse(data, content_type='application/octet-stream')
    response['Content-Disposition'] = f'inline; filename="{original_name}"'
    response['X-Original-Filename'] = original_name
    return response


# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

@require_POST
def cleanup(request):
    """Clean up the current session's work directory."""
    work_dir = request.session.get('dt_work_dir')
    if work_dir and os.path.isdir(work_dir):
        import shutil
        shutil.rmtree(work_dir, ignore_errors=True)

    # Clear session keys
    for key in ['dt_work_dir', 'dt_session_id', 'dt_blast_task_id', 'dt_uniprot_task_id']:
        request.session.pop(key, None)

    # Also clean old work directories
    clear_temp_directory(settings.WORK_DIRECTORY)

    return JsonResponse({'success': True})
