"""
Django views for the Heatmap Visualization web app.
"""
from __future__ import annotations

import json
import os
import uuid
import base64

from django.http import JsonResponse, HttpResponse
from django.shortcuts import render
from django.views.decorators.http import require_POST, require_GET, require_http_methods

from utils.lazy_import import LazyModule

# pandas and data_processor (which pulls in numpy) are heavy; the URLconf
# import chain loads this module on the very first request handled after
# boot (probe or otherwise), so both are deferred until a view actually
# uses them instead of paying that cost before the app can serve anything.
pd = LazyModule('pandas')
data_processor = LazyModule('heatmap_viz.services.data_processor')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_work_dir(request) -> str:
    work_dir = request.session.get('hm_work_dir')
    if work_dir and os.path.isdir(work_dir):
        return work_dir
    work_dir = os.path.join('/tmp', 'hm_' + str(uuid.uuid4()))
    os.makedirs(work_dir, mode=0o777, exist_ok=True)
    request.session['hm_work_dir'] = work_dir
    return work_dir


def _save_df(work_dir: str, name: str, df: pd.DataFrame):
    df.to_pickle(os.path.join(work_dir, f'{name}.pkl'))


def _load_df(work_dir: str, name: str) -> pd.DataFrame | None:
    path = os.path.join(work_dir, f'{name}.pkl')
    return pd.read_pickle(path) if os.path.exists(path) else None


def _save_json(work_dir: str, name: str, data):
    with open(os.path.join(work_dir, f'{name}.json'), 'w') as f:
        json.dump(data, f)


def _load_json(work_dir: str, name: str):
    path = os.path.join(work_dir, f'{name}.json')
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return None


# ---------------------------------------------------------------------------
# Main page
# ---------------------------------------------------------------------------

def index(request):
    return render(request, 'peptide/heatmap.html')


# ---------------------------------------------------------------------------
# Upload endpoint
# ---------------------------------------------------------------------------

@require_POST
def upload(request):
    """
    Accept merged CSV + optional FASTA file.
    POST body: multipart with 'merged_file' and optionally 'fasta_file'.
    """
    merged_file = request.FILES.get('merged_file')
    if not merged_file:
        return JsonResponse({'error': 'No merged data file uploaded.'}, status=400)

    df, group_data_dict, protein_dict, col_order, err = data_processor.load_merged_file(
        merged_file, merged_file.name
    )
    if err:
        return JsonResponse({'error': f'Could not load file: {err}'}, status=400)

    # Optional FASTA file – update protein_dict with sequences
    fasta_file = request.FILES.get('fasta_file')
    fasta_msg = None
    if fasta_file:
        merged_pids = set(df['Protein'].unique()) if 'Protein' in df.columns else None
        fasta_proteins, fasta_err = data_processor.load_fasta_file(fasta_file, merged_pids)
        if fasta_err:
            fasta_msg = f'FASTA warning: {fasta_err}'
        else:
            for pid, pdata in fasta_proteins.items():
                if pid in protein_dict:
                    protein_dict[pid].update(pdata)
                else:
                    protein_dict[pid] = pdata

    work_dir = _get_work_dir(request)
    _save_df(work_dir, 'merged_df', df)
    _save_json(work_dir, 'group_data_dict', group_data_dict)
    _save_json(work_dir, 'protein_dict', protein_dict)
    _save_json(work_dir, 'col_order', col_order)

    options = data_processor.get_selector_options(df, group_data_dict, protein_dict, col_order)

    warnings = []
    if fasta_msg:
        warnings.append(fasta_msg)

    return JsonResponse({
        'success': True,
        'filename': merged_file.name,
        'rows': len(df),
        'proteins': options['proteins'],
        'var_keys': options['var_keys'],
        'has_functions': options['has_functions'],
        'var_replicates': options['var_replicates'],
        'has_replicates': options['has_replicates'],
        'warnings': warnings,
    })


# ---------------------------------------------------------------------------
# Transfer from Data Transformation
# ---------------------------------------------------------------------------

@require_POST
def transfer_from_dt(request):
    """Load merged_df directly from the data transformation session."""
    import io
    dt_work_dir = request.session.get('dt_work_dir')
    if not dt_work_dir or not os.path.isdir(dt_work_dir):
        return JsonResponse(
            {'error': 'No data transformation session found. Please complete data transformation first.'},
            status=400,
        )
    dt_pkl = os.path.join(dt_work_dir, 'merged_df.pkl')
    if not os.path.exists(dt_pkl):
        return JsonResponse(
            {'error': 'No processed data found. Please click "Process Data" in the wizard first.'},
            status=400,
        )

    df = pd.read_pickle(dt_pkl)
    buf = io.BytesIO()
    df.to_csv(buf, index=False)
    buf.seek(0)

    df_hm, group_data_dict, protein_dict, col_order, err = data_processor.load_merged_file(buf, 'merged_df.csv')
    if err:
        return JsonResponse({'error': f'Could not process data: {err}'}, status=400)

    # protein_dict above is rebuilt from the merged CSV, which doesn't carry a
    # signal-peptide column — pull signal_end straight from Data Transformation's
    # own protein_dict.json (populated by its UniProt fetch) so "Strip start
    # sequence" has an auto-detected default after a transfer.
    dt_protein_dict = _load_json(dt_work_dir, 'protein_dict') or {}
    for pid, dt_info in dt_protein_dict.items():
        if pid in protein_dict and dt_info.get('signal_end') is not None:
            protein_dict[pid]['signal_end'] = dt_info['signal_end']

    # Carry over an optional FASTA provided during Data Transformation so the
    # heatmap has protein sequences loaded by default (mirrors a manual FASTA
    # upload on this page). Rebuilding protein_dict from the merged CSV above
    # means that when Data Transformation is reset — the wizard's "Start Over"
    # removes its whole work directory, including fasta_file.orig — this block
    # is skipped and any previously carried FASTA is cleared on the next transfer.
    fasta_filename = None
    fasta_warnings = []
    dt_fasta_path = os.path.join(dt_work_dir, 'fasta_file.orig')
    if os.path.exists(dt_fasta_path):
        merged_pids = set(df_hm['Protein'].unique()) if 'Protein' in df_hm.columns else None
        try:
            with open(dt_fasta_path, 'rb') as fh:
                fasta_proteins, fasta_err = data_processor.load_fasta_file(fh, merged_pids)
        except Exception as e:  # pragma: no cover - defensive
            fasta_proteins, fasta_err = {}, str(e)
        if fasta_err:
            fasta_warnings.append(f'FASTA warning: {fasta_err}')
        else:
            for pid, pdata in fasta_proteins.items():
                if pid in protein_dict:
                    protein_dict[pid].update(pdata)
                else:
                    protein_dict[pid] = pdata
        dt_names = _load_json(dt_work_dir, 'uploaded_file_names') or {}
        fasta_filename = dt_names.get('fasta_file') or 'Data Transformation FASTA'

    work_dir = _get_work_dir(request)
    _save_df(work_dir, 'merged_df', df_hm)
    _save_json(work_dir, 'group_data_dict', group_data_dict)
    _save_json(work_dir, 'protein_dict', protein_dict)
    _save_json(work_dir, 'col_order', col_order)

    options = data_processor.get_selector_options(df_hm, group_data_dict, protein_dict, col_order)
    return JsonResponse({
        'success': True,
        'filename': 'Data Transformation output',
        'rows': len(df_hm),
        'proteins': options['proteins'], 'var_keys': options['var_keys'],
        'has_functions': options['has_functions'],
        'var_replicates': options['var_replicates'],
        'has_replicates': options['has_replicates'],
        'fasta_filename': fasta_filename,
        'warnings': fasta_warnings,
    })


# ---------------------------------------------------------------------------
# Transfer from Data Analysis
# ---------------------------------------------------------------------------

@require_POST
def transfer_from_da(request):
    """Load merged_df directly from the data analysis session."""
    import io
    da_work_dir = request.session.get('da_work_dir')
    if not da_work_dir or not os.path.isdir(da_work_dir):
        return JsonResponse(
            {'error': 'No data analysis session found. Please upload a file in Data Analysis first.'},
            status=400,
        )
    da_pkl = os.path.join(da_work_dir, 'merged_df.pkl')
    if not os.path.exists(da_pkl):
        return JsonResponse(
            {'error': 'No processed data found in Data Analysis session.'},
            status=400,
        )

    df = pd.read_pickle(da_pkl)
    buf = io.BytesIO()
    df.to_csv(buf, index=False)
    buf.seek(0)

    df_hm, group_data_dict, protein_dict, col_order, err = data_processor.load_merged_file(buf, 'merged_df.csv')
    if err:
        return JsonResponse({'error': f'Could not process data: {err}'}, status=400)

    work_dir = _get_work_dir(request)
    _save_df(work_dir, 'merged_df', df_hm)
    _save_json(work_dir, 'group_data_dict', group_data_dict)
    _save_json(work_dir, 'protein_dict', protein_dict)
    _save_json(work_dir, 'col_order', col_order)

    options = data_processor.get_selector_options(df_hm, group_data_dict, protein_dict, col_order)
    return JsonResponse({
        'success': True,
        'filename': 'Data Analysis output',
        'rows': len(df_hm),
        'proteins': options['proteins'], 'var_keys': options['var_keys'],
        'has_functions': options['has_functions'],
        'var_replicates': options['var_replicates'],
        'has_replicates': options['has_replicates'], 'warnings': [],
    })


# ---------------------------------------------------------------------------
# Fetch protein sequence (UniProt)
# ---------------------------------------------------------------------------

@require_POST
def fetch_sequence(request):
    """
    Fetch sequence for a protein from UniProt if not already available.
    POST body (JSON): {"protein_id": "P12345"}
    """
    try:
        body = json.loads(request.body)
        protein_id = body.get('protein_id', '').strip()
    except Exception:
        return JsonResponse({'error': 'Invalid JSON body.'}, status=400)

    if not protein_id:
        return JsonResponse({'error': 'protein_id required.'}, status=400)

    work_dir = request.session.get('hm_work_dir')
    if not work_dir:
        return JsonResponse({'error': 'No session. Please upload a file first.'}, status=400)

    protein_dict = _load_json(work_dir, 'protein_dict') or {}

    # Already have sequence?
    existing = protein_dict.get(protein_id, {})
    if existing.get('sequence'):
        return JsonResponse({'success': True, 'cached': True, 'protein_id': protein_id})

    seq, signal_end = data_processor.fetch_sequence_from_uniprot(
        protein_id, species=existing.get('species'))
    if not seq:
        return JsonResponse({'error': f'Could not fetch sequence for {protein_id}.'}, status=404)

    entry = protein_dict.setdefault(protein_id, {})
    entry['sequence'] = seq
    entry['signal_end'] = signal_end
    _save_json(work_dir, 'protein_dict', protein_dict)

    return JsonResponse({
        'success': True, 'cached': False,
        'protein_id': protein_id, 'length': len(seq),
    })


# ---------------------------------------------------------------------------
# Apply a FASTA to already-loaded session data
# ---------------------------------------------------------------------------

@require_POST
def apply_fasta(request):
    """Apply an uploaded FASTA to the merged data already in the session.

    The /upload/ endpoint attaches a FASTA only when the merged file is POSTed
    with it. But when the merged data was loaded without a client-side File
    object — e.g. transferred from Data Transformation — re-POSTing the merged
    file isn't possible, so a FASTA picked afterwards would otherwise be dropped
    and the app would fall back to (often failing) UniProt lookups. This merges
    the FASTA sequences into the saved protein_dict and returns refreshed
    protein options so the selector's has_sequence flags update in place.
    """
    fasta_file = request.FILES.get('fasta_file')
    if not fasta_file:
        return JsonResponse({'error': 'No FASTA file uploaded.'}, status=400)

    work_dir = request.session.get('hm_work_dir')
    if not work_dir:
        return JsonResponse({'error': 'No session. Please load data first.'}, status=400)

    merged_df = _load_df(work_dir, 'merged_df')
    if merged_df is None or 'Protein' not in merged_df.columns:
        return JsonResponse({'error': 'No merged data in session. Please load data first.'}, status=400)

    protein_dict = _load_json(work_dir, 'protein_dict') or {}
    group_data_dict = _load_json(work_dir, 'group_data_dict') or {}
    col_order = _load_json(work_dir, 'col_order') or []

    merged_pids = set(merged_df['Protein'].unique())
    fasta_proteins, fasta_err = data_processor.load_fasta_file(fasta_file, merged_pids)
    if fasta_err:
        return JsonResponse({'error': f'FASTA error: {fasta_err}'}, status=400)

    applied = 0
    for pid, pdata in fasta_proteins.items():
        if pid in protein_dict:
            protein_dict[pid].update(pdata)
        else:
            protein_dict[pid] = pdata
        if pdata.get('sequence'):
            applied += 1

    _save_json(work_dir, 'protein_dict', protein_dict)

    options = data_processor.get_selector_options(
        merged_df, group_data_dict, protein_dict, col_order
    )
    return JsonResponse({
        'success': True,
        'filename': fasta_file.name,
        'sequences_applied': applied,
        'proteins': options['proteins'],
    })


# ---------------------------------------------------------------------------
# Get specific options (for bio_or_pep dropdown)
# ---------------------------------------------------------------------------

@require_POST
def get_specific_options(request):
    """
    Return options for specific_select_multiple based on bio_or_pep value.
    POST body (JSON): {"bio_or_pep": "2", "selected_proteins": [...]}
    """
    try:
        body = json.loads(request.body)
        bio_or_pep = body.get('bio_or_pep', 'no')
        selected_proteins = body.get('selected_proteins', [])
        selected_var_keys = body.get('selected_var_keys', [])
    except Exception:
        return JsonResponse({'error': 'Invalid JSON body.'}, status=400)

    work_dir = request.session.get('hm_work_dir')
    if not work_dir:
        return JsonResponse({'error': 'No session. Please upload a file first.'}, status=400)

    merged_df = _load_df(work_dir, 'merged_df')
    if merged_df is None:
        return JsonResponse({'error': 'Session data missing.'}, status=400)

    options = data_processor.get_specific_options(
        merged_df, bio_or_pep,
        selected_proteins or None,
        selected_var_keys or None,
    )
    return JsonResponse({'options': options})


# ---------------------------------------------------------------------------
# "Strip start sequence" — UniProt signal-peptide suggestion
# ---------------------------------------------------------------------------

@require_POST
def signal_peptide_info(request):
    """
    Return the UniProt-annotated signal-peptide suggestion (position + amino
    acid letters) for each requested protein, so the "Strip start sequence"
    UI can show the user what it would auto-detect next to the manual
    residue-count override.
    POST body (JSON): {"selected_proteins": ["P12345", ...]}
    """
    try:
        body = json.loads(request.body)
        selected_proteins = body.get('selected_proteins', [])
    except Exception:
        return JsonResponse({'error': 'Invalid JSON body.'}, status=400)

    work_dir = request.session.get('hm_work_dir')
    if not work_dir:
        return JsonResponse({'error': 'No session. Please upload a file first.'}, status=400)

    protein_dict = _load_json(work_dir, 'protein_dict') or {}
    suggestions = data_processor.get_signal_peptide_suggestions(protein_dict, selected_proteins)
    return JsonResponse({'suggestions': suggestions})


# ---------------------------------------------------------------------------
# Plot endpoint
# ---------------------------------------------------------------------------

@require_POST
def plot(request):
    """
    Generate heatmap images.
    POST body (JSON): {selected_proteins, selected_var_keys, plot_params}
    Returns base64 PNG images.
    """
    try:
        body = json.loads(request.body)
    except Exception:
        return JsonResponse({'error': 'Invalid JSON body.'}, status=400)

    work_dir = request.session.get('hm_work_dir')
    if not work_dir:
        return JsonResponse({'error': 'No session. Please upload a file first.'}, status=400)

    merged_df = _load_df(work_dir, 'merged_df')
    group_data_dict = _load_json(work_dir, 'group_data_dict')
    protein_dict = _load_json(work_dir, 'protein_dict')

    if merged_df is None or group_data_dict is None:
        return JsonResponse({'error': 'Session data missing. Please re-upload your file.'}, status=400)

    # Belt-and-suspenders check: load_merged_file() already rejects a dataset
    # with no start/end columns at all, but a dataset that reached this session
    # via a transfer path could in principle carry 'start'/'end' columns that
    # are present but entirely empty. Either way, a sequence heatmap cannot be
    # generated without real peptide position data, so refuse outright rather
    # than attempt a plot and fail confusingly deeper in the pipeline.
    if (
        'start' not in merged_df.columns or 'end' not in merged_df.columns
        or merged_df['start'].isna().all() or merged_df['end'].isna().all()
    ):
        return JsonResponse({
            'error': 'Missing peptide position data.',
            'message': (
                'This dataset has no usable start/end (peptide position) values, '
                'so a sequence heatmap cannot be generated. Upload data that '
                'includes peptide start and end positions.'
            ),
        }, status=422)

    selected_proteins = body.get('selected_proteins', [])
    selected_var_keys = body.get('selected_var_keys', [])
    plot_params = body.get('plot_params', {})

    if not selected_proteins or not selected_var_keys:
        return JsonResponse({'error': 'Please select at least one protein and one variable.'}, status=400)

    # Check sequences – if missing, report which need fetching
    missing_seqs = [
        pid for pid in selected_proteins
        if not protein_dict.get(pid, {}).get('sequence', '')
    ]
    if missing_seqs:
        return JsonResponse({
            'error': 'Missing protein sequences.',
            'missing_sequences': missing_seqs,
            'message': (
                f'Sequences for {", ".join(missing_seqs)} are not available. '
                'Please upload a FASTA file or enable UniProt lookup.'
            ),
        }, status=422)

    # Build available_data_variables_dict
    available, build_msgs = data_processor.build_available_data_variables(
        merged_df, protein_dict, group_data_dict,
        selected_proteins, selected_var_keys,
        strip_start_sequence=plot_params.get('strip_start_sequence', False),
        strip_start_manual=plot_params.get('strip_start_manual'),
    )

    if not available:
        return JsonResponse({
            'error': 'Could not build plot data.',
            'messages': build_msgs,
        }, status=400)

    # Generate heatmap
    portrait_b64, landscape_b64, compact_json, \
        portrait_interactive, landscape_interactive, plot_msgs = data_processor.generate_heatmap(
            available, plot_params
        )

    all_msgs = build_msgs + plot_msgs

    if portrait_b64 is None and landscape_b64 is None and compact_json is None \
            and portrait_interactive is None and landscape_interactive is None:
        return JsonResponse({
            'error': 'No heatmap generated.',
            'messages': all_msgs,
        }, status=400)

    return JsonResponse({
        'success': True,
        'portrait':              portrait_b64,
        'landscape':             landscape_b64,
        'compact':               compact_json,
        'portrait_interactive':  portrait_interactive,
        'landscape_interactive': landscape_interactive,
        'messages': all_msgs,
    })


# ---------------------------------------------------------------------------
# Save / download compact interactive plot as HTML
# ---------------------------------------------------------------------------

@require_POST
def save_compact_plot(request):
    """Save a client-generated compact plot HTML to the session work dir."""
    try:
        data = json.loads(request.body)
        html_content = data.get('html', '')
    except Exception:
        return JsonResponse({'error': 'Invalid body.'}, status=400)

    if not html_content:
        return JsonResponse({'error': 'No HTML content provided.'}, status=400)

    work_dir = _get_work_dir(request)
    plot_path = os.path.join(work_dir, 'compact_plot.html')
    with open(plot_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    return JsonResponse({'success': True})


@require_GET
def download_compact_plot(request):
    """Return the previously saved compact plot as a self-contained HTML file."""
    work_dir = request.session.get('hm_work_dir')
    if not work_dir:
        return HttpResponse('No session data found.', status=400)

    plot_path = os.path.join(work_dir, 'compact_plot.html')
    if not os.path.exists(plot_path):
        return HttpResponse('No compact plot has been generated yet.', status=404)

    with open(plot_path, 'rb') as f:
        response = HttpResponse(f.read(), content_type='text/html')
        response['Content-Disposition'] = 'attachment; filename="heatmap_compact.html"'
        return response


# ---------------------------------------------------------------------------
# Download plot as PNG
# ---------------------------------------------------------------------------

@require_POST
def download_plot(request):
    """
    Accept base64 PNG from frontend, decode, return as file attachment.
    POST body (JSON): {"image_b64": "...", "orientation": "portrait"}
    """
    try:
        body = json.loads(request.body)
        b64 = body.get('image_b64', '')
        orientation = body.get('orientation', 'portrait')
    except Exception:
        return HttpResponse('Invalid request.', status=400)

    if not b64:
        return HttpResponse('No image data.', status=400)

    img_bytes = base64.b64decode(b64)
    filename = f'heatmap_{orientation}.png'
    response = HttpResponse(img_bytes, content_type='image/png')
    response['Content-Disposition'] = f'attachment; filename="{filename}"'
    return response


@require_POST
def download_static_image(request):
    """Server-side Kaleido export of the currently shown interactive heatmap.

    POST JSON: {figure_json, format: 'png'|'svg', filename, width?, height?}.
    PNG is rendered at publication DPI (see utils.plotly_export.PUBLICATION_DPI);
    SVG is true vector. Returns the file as an attachment.
    """
    from utils import plotly_export

    try:
        body = json.loads(request.body)
        figure_json = body.get('figure_json')
        fmt = (body.get('format') or 'png').lower()
        raw_name = body.get('filename') or 'heatmap'
        width = body.get('width') or None
        height = body.get('height') or None
    except Exception:
        return HttpResponse('Invalid request.', status=400)

    if not figure_json:
        return HttpResponse('No figure data.', status=400)

    # Sanitize the download filename (header-safe).
    filename = ''.join(c if c.isalnum() or c in ('-', '_', '.') else '_' for c in str(raw_name)) or 'heatmap'

    # Below-axis legends are positioned in paper fractions from a pixel-width
    # estimate, so re-lay them for the width this export actually renders at
    # (no-op for right-legend figures — see repack_below_legends).
    from heatmap_viz.services.heatmap_renderer import repack_below_legends
    figure_json = repack_below_legends(figure_json, width)

    try:
        if fmt == 'svg':
            data = plotly_export.svg_bytes(figure_json, width=width, height=height)
            content_type, ext = 'image/svg+xml', 'svg'
        else:
            data = plotly_export.png_bytes(figure_json, width=width, height=height)
            content_type, ext = 'image/png', 'png'
    except Exception as exc:
        return HttpResponse(f'Static export failed: {exc}', status=500)

    response = HttpResponse(data, content_type=content_type)
    response['Content-Disposition'] = f'attachment; filename="{filename}.{ext}"'
    return response
