"""
Django views for the Data Analysis web app.
"""
from __future__ import annotations

import json
import os
import uuid

from django.http import JsonResponse, HttpResponse
from django.shortcuts import render
from django.views.decorators.http import require_POST, require_GET

from utils.lazy_import import LazyModule

# pandas/numpy/scipy are heavy; the URLconf import chain loads this module
# on the very first request handled after boot (probe or otherwise), so
# they're deferred until a view actually uses them.
pd = LazyModule('pandas')
data_processor = LazyModule('data_analysis.services.data_processor')
plotter = LazyModule('data_analysis.services.plotter')


# ---------------------------------------------------------------------------
# Helpers – session-scoped work directory
# ---------------------------------------------------------------------------

def _get_work_dir(request) -> str:
    work_dir = request.session.get('da_work_dir')
    if work_dir and os.path.isdir(work_dir):
        return work_dir
    work_dir = os.path.join('/tmp', 'da_' + str(uuid.uuid4()))
    os.makedirs(work_dir, mode=0o777, exist_ok=True)
    request.session['da_work_dir'] = work_dir
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
    return render(request, 'peptide/data_analysis.html')


# ---------------------------------------------------------------------------
# Upload endpoint
# ---------------------------------------------------------------------------

@require_POST
def upload(request):
    """
    Accept a merged data CSV file, process it, return selector options.
    POST body: multipart with 'merged_file'
    """
    merged_file = request.FILES.get('merged_file')
    if not merged_file:
        return JsonResponse({'error': 'No file uploaded.'}, status=400)

    df, err = data_processor.load_file(merged_file, merged_file.name)
    if err or df is None:
        return JsonResponse({'error': f'Could not load file: {err}'}, status=400)

    df, col_warnings = data_processor.validate_columns(df)

    group_data_dict, df, group_warning = data_processor.process_group_data(df)
    if not group_data_dict:
        return JsonResponse(
            {'error': 'No group columns found. Expected Avg_* columns or columns with Grouped: format.'},
            status=400,
        )

    protein_dict = data_processor.extract_protein_dict(df)

    # Persist to session work dir
    work_dir = _get_work_dir(request)
    _save_df(work_dir, 'merged_df', df)
    _save_json(work_dir, 'group_data_dict', group_data_dict)
    _save_json(work_dir, 'protein_dict', protein_dict)

    options = data_processor.get_selector_options(df, group_data_dict, protein_dict)
    warnings = col_warnings + ([group_warning] if group_warning else [])

    return JsonResponse({
        'success': True,
        'filename': merged_file.name,
        'rows': len(df),
        'columns': len(df.columns),
        'groups': options['groups'],
        'proteins': options['proteins'],
        'protein_ids': options['protein_ids'],
        'functions': options['functions'],
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
    df, col_warnings = data_processor.validate_columns(df)
    group_data_dict, df, group_warning = data_processor.process_group_data(df)
    if not group_data_dict:
        return JsonResponse({'error': 'No group columns found in transformed data.'}, status=400)
    protein_dict = data_processor.extract_protein_dict(df)

    work_dir = _get_work_dir(request)
    _save_df(work_dir, 'merged_df', df)
    _save_json(work_dir, 'group_data_dict', group_data_dict)
    _save_json(work_dir, 'protein_dict', protein_dict)

    options = data_processor.get_selector_options(df, group_data_dict, protein_dict)
    warnings = col_warnings + ([group_warning] if group_warning else [])
    return JsonResponse({
        'success': True,
        'filename': 'Data Transformation output',
        'rows': len(df), 'columns': len(df.columns),
        'groups': options['groups'], 'proteins': options['proteins'],
        'protein_ids': options['protein_ids'], 'functions': options['functions'],
        'has_functions': options['has_functions'],
        'var_replicates': options['var_replicates'],
        'has_replicates': options['has_replicates'], 'warnings': warnings,
    })


# ---------------------------------------------------------------------------
# Plot endpoint
# ---------------------------------------------------------------------------

@require_POST
def plot(request):
    """
    Generate a Plotly chart and return JSON spec.
    POST body (JSON): plot parameters
    """
    try:
        params = json.loads(request.body)
    except Exception:
        return JsonResponse({'error': 'Invalid JSON body.'}, status=400)

    work_dir = request.session.get('da_work_dir')
    if not work_dir:
        return JsonResponse({'error': 'No data uploaded. Please upload a file first.'}, status=400)

    merged_df = _load_df(work_dir, 'merged_df')
    group_data_dict = _load_json(work_dir, 'group_data_dict')
    protein_dict = _load_json(work_dir, 'protein_dict')

    if merged_df is None or group_data_dict is None:
        return JsonResponse({'error': 'Session data missing. Please re-upload your file.'}, status=400)

    fig_json, warnings = plotter.generate_plot(merged_df, group_data_dict, protein_dict, params)

    if fig_json is None:
        return JsonResponse({'error': 'Failed to generate plot.', 'warnings': warnings}, status=400)

    return JsonResponse({'figure': json.loads(fig_json), 'warnings': warnings})


# ---------------------------------------------------------------------------
# Download plot as HTML
# ---------------------------------------------------------------------------

def _sanitize_filename(raw_name, fallback):
    """Same title-derived naming scheme used across all plot export formats."""
    cleaned = ''.join(c if c.isalnum() or c in ('-', '_', '.') else '_' for c in str(raw_name))
    return cleaned or fallback


@require_GET
def download_plot(request):
    """Return a previously generated plot as a self-contained HTML file."""
    work_dir = request.session.get('da_work_dir')
    if not work_dir:
        return HttpResponse('No session data found.', status=400)

    plot_path = os.path.join(work_dir, 'current_plot.html')
    if not os.path.exists(plot_path):
        return HttpResponse('No plot has been generated yet.', status=404)

    filename = _sanitize_filename(request.GET.get('filename') or 'data_analysis_plot', 'data_analysis_plot')

    with open(plot_path, 'rb') as f:
        response = HttpResponse(f.read(), content_type='text/html')
        response['Content-Disposition'] = f'attachment; filename="{filename}.html"'
        return response


@require_POST
def download_static_image(request):
    """Server-side Kaleido export of the current Data Analysis figure.

    POST JSON: {figure_json, format: 'png'|'svg', filename, width?, height?}.
    PNG is rendered at publication DPI (see utils.plotly_export.PUBLICATION_DPI);
    SVG is true vector. Returns the file as an attachment.
    """
    from utils import plotly_export

    try:
        body = json.loads(request.body)
        figure_json = body.get('figure_json')
        fmt = (body.get('format') or 'png').lower()
        raw_name = body.get('filename') or 'peptiline_plot'
        width = body.get('width') or None
        height = body.get('height') or None
    except Exception:
        return HttpResponse('Invalid request.', status=400)

    if not figure_json:
        return HttpResponse('No figure data.', status=400)

    filename = _sanitize_filename(raw_name, 'peptiline_plot')

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


@require_POST
def save_plot(request):
    """Save the current plot to disk so it can be downloaded."""
    try:
        data = json.loads(request.body)
        html_content = data.get('html', '')
    except Exception:
        return JsonResponse({'error': 'Invalid body.'}, status=400)

    work_dir = _get_work_dir(request)
    plot_path = os.path.join(work_dir, 'current_plot.html')
    with open(plot_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    return JsonResponse({'success': True})
