"""
Heatmap rendering functions for heatmap_visualization.ipynb.

Importing this module does NOT trigger matplotlib or seaborn loading.
Those heavy libraries are imported lazily inside the rendering functions
that actually use them — on the first plot call, not on notebook startup.
"""

import html
import json
import logging
import pandas as pd
import numpy as np
import re
import copy
import warnings
from itertools import cycle

logger = logging.getLogger(__name__)


def _log_validation(message: str, level: str = "warning") -> None:
    """Log validation/warning messages (replaces legacy IPython display)."""
    if level == "error":
        logger.error(message)
    elif level == "info":
        logger.info(message)
    else:
        logger.warning(message)


def numeric_column(frame: pd.DataFrame, col) -> pd.Series:
    """
    Return ``frame[col]`` coerced to numeric, always as a 1-d Series.

    Peptide-interval column labels ("147-158 SANDGRGDSVAY") are not unique once
    more than a couple of proteins are selected, because the same peptide can map
    to several parent proteins. When a label is duplicated, ``frame[col]`` returns
    a DataFrame rather than a Series, and ``pd.to_numeric`` rejects it with
    "arg must be a list, tuple, 1-d array, or Series".

    Duplicates describe the same peptide at the same interval, so they are
    collapsed with a row-wise max: that preserves both peak abundance and
    non-zero detection, which is all any caller uses these values for.
    """
    selected = frame[col]
    if isinstance(selected, pd.DataFrame):
        selected = selected.apply(pd.to_numeric, errors='coerce').max(axis=1)
    return pd.to_numeric(selected, errors='coerce')

# ---------------------------------------------------------------------------
# Rendering constants (formerly imported from notebooks/_settings.py)
# ---------------------------------------------------------------------------
plot_heatmap     = 'yes'
plot_zero        = 'no'
chuck_size       = 78
legend_title     = [
    'Sample Type:',
    'Peptide Counts:',
    'Bioactivity Function:',
    'Peptide Interval:',
    'Average Abundance:',
]
# ---------------------------------------------------------------------------
# State globals — set by process_available_data()/update_plot() at runtime,
# read by the visualize_* functions in this same module.
# ---------------------------------------------------------------------------
y_ticks          = []
max_sequence_length = 0
max_count        = 0
num_unique_count = 0
num_colors       = 0

# Typeface for the interactive (Plotly) heatmap figures. Kept identical to the
# Data Analysis dashboard (data_analysis/services/plotter.py FONT_FAMILY) so the
# two dashboards share one look. Applied as layout.font so it cascades to title,
# axes and legend.
PLOTLY_FONT_FAMILY = 'Arial, Helvetica, sans-serif'

# ---------------------------------------------------------------------------
# Legend placement
# ---------------------------------------------------------------------------
# 'right'  — the historical layout: ONE combined Plotly legend to the right of
#            the plot, every group (Sample Type, Peptide Counts, Average
#            Abundance / Standardized Mean Difference / log2 Fold Change)
#            stacked inside it, text reading downwards.
# 'below'  — each group becomes its OWN boxed legend unit placed under the
#            x-axis (i.e. under the protein-sequence label), one per row,
#            centred, stacked above one another. Nothing is reserved to the
#            right, so the heatmap keeps the full figure width. (A side-by-side
#            packing was tried and dropped — it read poorly once wrapped onto
#            multiple rows, so 'below' is stacked-only now.)
LEGEND_POSITIONS = ('right', 'below')

# Geometry of the 'below' legend strip, in pixels, tuned at LEGEND_FONT_REF_PX
# and scaled from there by the actual font_size_legend (see _below_legend_layout)
# so a much larger legend font gets a proportionally taller box instead of
# overflowing it. Plotly anchors each legend box for you but never lays
# MULTIPLE legends out relative to each other, so the units are positioned here
# from an estimate of how wide each one draws. The estimate is in pixels while
# Plotly positions in paper fractions, so the layout needs the width of the
# plotting area to convert between them:
#   * the browser knows it exactly — heatmap.html re-packs on render and resize
#     (repackBelowLegends), and the static-export view re-packs server-side with
#     the width the client sends (see repack_below_legends below);
#   * LEGEND_REF_WIDTH_PX is only the fallback used to build the first layout,
#     before any of that runs.
LEGEND_REF_WIDTH_PX = 900     # fallback width when the real one isn't known yet
LEGEND_FONT_REF_PX = 14       # font size the pixel constants below are tuned for
LEGEND_UNIT_GAP_PX = 30       # gap between neighbouring units in a row
LEGEND_ROW_PX = 40            # vertical room one row of units needs (single inline
                               # line now that the title sits left of its entries
                               # instead of stacked above them) at LEGEND_FONT_REF_PX
LEGEND_AXIS_ZONE_PX = 90      # x-axis ticks + title, cleared before the strip, at
                               # LEGEND_FONT_REF_PX — generous enough that a
                               # wrapped/large x-axis title never collides with it
LEGEND_BASE_BOTTOM_PX = 80    # bottom margin below the whole strip
LEGEND_OUTLINE_PAD_PX = 4     # padding between the strip's content and its
                               # single outline, at LEGEND_FONT_REF_PX (see
                               # _below_legend_geometry)


def normalize_legend_position(value) -> str:
    """Coerce a user-supplied legend position to one of :data:`LEGEND_POSITIONS`."""
    text = str(value or '').strip().lower().replace('_', '-').replace(' ', '-')
    if text in ('below', 'below-stacked', 'bottom', 'bottom-stacked', 'stacked', 'under'):
        return 'below'
    return 'right'


def _estimate_legend_unit_px(title, labels, font_size_legend) -> float:
    """Approximate drawn width of one horizontal legend unit, in pixels.

    The title sits inline to the LEFT of its entries (``title.side='left'``),
    not stacked above them, so the unit's width is the title's width PLUS the
    entries' width rather than the max of the two.
    """
    char_px = 0.58 * font_size_legend
    items_px = sum(40 + len(str(lbl)) * char_px for lbl in labels)
    title_px = len(str(title)) * char_px + 20
    return title_px + items_px + 24            # + box padding / border


def _make_legend_placer(legend_position, font_size_legend):
    """Build the (kwargs-factory, registry) pair used to place legend entries.

    Returns ``(legend_kwargs, units)`` where ``legend_kwargs(group_title,
    legendgroup, entry=None)`` yields the trace kwargs for one legend entry and
    ``units`` is the ordered ``{group title: {'ref', 'labels'}}`` registry —
    empty in 'right' mode, one entry per group (in first-drawn order, i.e. left
    to right) in the 'below' modes. ``entry`` is the visible item label, recorded
    so the unit's drawn width can be estimated. See :data:`LEGEND_POSITIONS`.
    """
    below = normalize_legend_position(legend_position) == 'below'
    units = {}

    def legend_kwargs(group_title, legendgroup, entry=None):
        if not below:
            return dict(
                legendgroup=legendgroup,
                legendgrouptitle=dict(text=group_title, font=dict(size=font_size_legend)),
            )
        unit = units.setdefault(
            group_title,
            {'ref': 'legend' if not units else f'legend{len(units) + 1}', 'labels': []},
        )
        if entry is not None:
            unit['labels'].append(entry)
        # The unit's own legend title carries the header, so the per-group title
        # is dropped — it would otherwise print twice in the same box. The
        # legendgroup goes too: inside a HORIZONTAL legend Plotly lays each
        # legendgroup out as its own vertical column, so keeping it would stack
        # the entries downwards and defeat the point of the below-axis strip.
        # Each unit is already one group, so nothing is lost.
        return dict(legend=unit['ref'])

    return legend_kwargs, units


def _below_legend_plan(units, font_size_legend):
    """Serialisable [{ref, title, width}] for every registered unit, in draw order."""
    return [
        {'ref': unit['ref'], 'title': title,
         'width': _estimate_legend_unit_px(title, unit['labels'], font_size_legend)}
        for title, unit in units.items()
    ]


def _pack_legend_rows(plan, avail_px, stacked=False):
    """Group the planned units into rows that each fit ``avail_px``.

    ``stacked`` forces one unit per row, which no width can crowd. Returns a
    list of rows, each a list of the plan's unit dicts.
    """
    rows, row, row_px = [], [], 0.0
    for unit in plan:
        too_wide = row_px + LEGEND_UNIT_GAP_PX + unit['width'] > avail_px
        if row and (stacked or too_wide):
            rows.append(row)
            row, row_px = [], 0.0
        row_px += unit['width'] + (LEGEND_UNIT_GAP_PX if row else 0)
        row.append(unit)
    if row:
        rows.append(row)
    return rows


def _outline_pad_px(plan, font_scale=1.0) -> float:
    """Padding between the strip's content and its single outline, in pixels.

    Three groups (Sample Type / Peptide Counts / Average Abundance-or-SMD-or-
    log2FC — the Portrait comparison case) pack noticeably tighter than one or
    two, so they get a little more breathing room. ``font_scale`` (see
    :func:`_below_legend_layout`) grows the padding along with the legend font.
    """
    base = LEGEND_OUTLINE_PAD_PX + (8 if len(plan) == 3 else 0)
    return base * font_scale


def _below_legend_geometry(plan, stacked, avail_w_px, top_px, base_height_px,
                           row_px=LEGEND_ROW_PX, axis_zone_px=LEGEND_AXIS_ZONE_PX,
                           outline_pad_px=None):
    """Where every 'below' legend unit goes, for a known plotting-area width.

    This is the single source of truth for the strip's arithmetic — the initial
    server-side layout, the browser's re-pack, and the static-export re-pack all
    call it (the browser through the JS mirror in heatmap.html). Returns
    ``(positions, height_px, bottom_margin_px, n_rows, bbox)`` where
    ``positions`` maps each legend ref to its ``(x, y)`` in Plotly paper
    fractions, and ``bbox`` is ``(x0, x1, y0, y1)`` — the single outline drawn
    around the WHOLE strip (all units, all rows) instead of a border per unit.

    ``row_px``, ``axis_zone_px`` and ``outline_pad_px`` default to the module
    constants (tuned at :data:`LEGEND_FONT_REF_PX`) but callers pass in
    font-scaled values — see :func:`_below_legend_layout`. Units are always one
    per row, centred: the 'below' position is stacked-only (a side-by-side
    packing was tried and dropped for looking crowded once wrapped).
    """
    rows = _pack_legend_rows(plan, avail_w_px, stacked)
    n_rows = max(len(rows), 1)
    strip_px = row_px * n_rows
    height = base_height_px + strip_px
    bottom = LEGEND_BASE_BOTTOM_PX + strip_px
    plot_h = max(height - top_px - bottom, 1)

    positions = {}
    x0, x1 = float('inf'), float('-inf')
    for r, row in enumerate(rows):
        span = sum(u['width'] for u in row) + LEGEND_UNIT_GAP_PX * (len(row) - 1)
        cursor = -span / 2      # centred cluster (one unit per row, so just that unit)
        y = -((axis_zone_px + r * row_px) / plot_h)
        for unit in row:
            positions[unit['ref']] = (0.5 + (cursor + unit['width'] / 2) / avail_w_px, y)
            x0 = min(x0, 0.5 + cursor / avail_w_px)
            x1 = max(x1, 0.5 + (cursor + unit['width']) / avail_w_px)
            cursor += unit['width'] + LEGEND_UNIT_GAP_PX

    pad = outline_pad_px if outline_pad_px is not None else _outline_pad_px(plan)
    pad_x = pad / avail_w_px
    pad_y = pad / plot_h
    bbox = (
        x0 - pad_x, x1 + pad_x,
        -(axis_zone_px / plot_h) + pad_y,                        # top
        -((axis_zone_px + n_rows * row_px) / plot_h) - pad_y,    # bottom
    )
    return positions, height, bottom, len(rows), bbox


def _below_legend_layout(plan, stacked, top_px, base_height_px, font_size_legend,
                         avail_w_px=None):
    """Full layout for the 'below' strip: legend objects, height, margin, meta.

    Returns ``(layout, height_px, bottom_margin_px, meta)``. ``meta`` is stored
    on ``layout.meta`` so the browser and the export path can re-pack the units
    once the real width is known — see :func:`repack_below_legends`.
    """
    # Row height, axis clearance and outline padding all grow with the legend
    # font — otherwise a large font_size_legend (e.g. 28 vs. the tuned-at 14)
    # overflows a box sized for the small default and crowds the x-axis title
    # above it.
    font_scale = font_size_legend / LEGEND_FONT_REF_PX
    row_px = LEGEND_ROW_PX * font_scale
    axis_zone_px = LEGEND_AXIS_ZONE_PX * font_scale
    outline_pad_px = _outline_pad_px(plan, font_scale)

    avail = avail_w_px or LEGEND_REF_WIDTH_PX
    positions, height, bottom, _rows, _bbox = _below_legend_geometry(
        plan, stacked, avail, top_px, base_height_px,
        row_px=row_px, axis_zone_px=axis_zone_px, outline_pad_px=outline_pad_px)

    layout = {}
    for unit in plan:
        x, y = positions[unit['ref']]
        layout[unit['ref']] = dict(
            # Title sits inline to the LEFT of its entries ("Peptide Counts: ⬤ 0  ⬤ 1…")
            # instead of stacked above them.
            title=dict(text=unit['title'],
                       font=dict(size=font_size_legend, color='black'),
                       side='left'),
            orientation='h',
            x=x, xanchor='center',
            y=y, yanchor='top',
            borderwidth=0,
            bgcolor='rgba(255,255,255,0.9)',
            font=dict(size=font_size_legend),
            itemsizing='constant',
        )

    meta = {'below_legend': {
        'units': plan,
        'stacked': bool(stacked),
        'gap_px': LEGEND_UNIT_GAP_PX,
        'row_px': row_px,
        'axis_zone_px': axis_zone_px,
        'base_bottom_px': LEGEND_BASE_BOTTOM_PX,
        'base_height_px': base_height_px,
        'top_px': top_px,
        'outline_pad_px': outline_pad_px,
    }}
    return layout, height, bottom, meta


def repack_below_legends(figure_json: str, width_px) -> str:
    """Re-lay the 'below' legend units for a figure being exported at ``width_px``.

    The server builds the first layout against :data:`LEGEND_REF_WIDTH_PX`
    because it cannot know the browser's width; a static export at a different
    width would otherwise place the units for the wrong figure — crowding them
    when it is narrower. Returns the figure JSON unchanged when it carries no
    'below' legend plan (every 'right'-legend figure, and every figure from the
    other dashboards) or when no width is supplied.
    """
    if not figure_json or not width_px:
        return figure_json
    try:
        fig = json.loads(figure_json) if isinstance(figure_json, str) else figure_json
        layout = fig.get('layout') or {}
        plan_meta = (layout.get('meta') or {}).get('below_legend')
        if not plan_meta:
            return figure_json

        margin = layout.get('margin') or {}
        avail = float(width_px) - float(margin.get('l', 0)) - float(margin.get('r', 0))
        if avail <= 0:
            return figure_json

        positions, height, bottom, _rows, _bbox = _below_legend_geometry(
            plan_meta['units'], plan_meta['stacked'], avail,
            plan_meta['top_px'], plan_meta['base_height_px'],
            row_px=plan_meta['row_px'], axis_zone_px=plan_meta['axis_zone_px'],
            outline_pad_px=plan_meta['outline_pad_px'])

        for ref, (x, y) in positions.items():
            legend = layout.setdefault(ref, {})
            legend['x'], legend['y'] = x, y
        layout['height'] = height
        margin['b'] = bottom
        layout['margin'] = margin
        fig['layout'] = layout
        return json.dumps(fig)
    except Exception as exc:      # never block an export over legend placement
        logger.warning('Could not re-pack below-axis legends: %s', exc)
        return figure_json


# Define functions outside of class
def update_filenames(input_filename_port, input_filename_land):
    # Append directory and .png only when necessary
    #updated_filename_port = f'{images_folder_name}/{input_filename_port}.png' if input_filename_port else None
    #updated_filename_land = f'{images_folder_name}/{input_filename_land}.png' if input_filename_land else None
    updated_filename_port = f'{input_filename_port}' if input_filename_port else None
    updated_filename_land = f'{input_filename_land}' if input_filename_land else None

    return updated_filename_land, updated_filename_port

def proceed_with_label_specific_options(selected_bio_or_pep, bio_or_pep):
    # Initialize selected_peptides and selected_functions to ensure they are always defined
    selected_peptides = []
    selected_functions = []

    # Check and handle different cases based on bio_or_pep value
    if bio_or_pep == '1':  # Assuming '1' indicates selection of peptides
        selected_peptides = list(selected_bio_or_pep) if isinstance(selected_bio_or_pep, (list, tuple)) else [
            selected_bio_or_pep]
    elif bio_or_pep == '2':  # Assuming '2' indicates selection of functions
        selected_functions = list(selected_bio_or_pep) if isinstance(selected_bio_or_pep, (list, tuple)) else [
            selected_bio_or_pep]
    # Optional: handle unexpected bio_or_pep values
    else:
        error_message = f"Unexpected value for bio_or_pep: {bio_or_pep}"
        _log_validation(error_message)
    #print(f" DEBUG: selected_peptides: {selected_peptides}")
    return selected_peptides, selected_functions

def get_interval_start(peptide_label):
    """
    Extract numeric start from labels like '193-207' or '193-207 (Oxidation)'.
    """
    base = peptide_label.split()[0]  # take '193-207' from '193-207 (Oxidation)'
    return int(base.split('-')[0])

"""-----------------Export_Function---------------------------------"""


def _log_base_fn(log_base):
    """Return the numpy log function for base 10 (default) or 2."""
    return np.log2 if int(log_base) == 2 else np.log10


def _log_sub(log_base):
    """Subscript for axis/legend labels: '10' or '2'."""
    return '2' if int(log_base) == 2 else '10'


def _log_value(v, log_base):
    """Base-aware scalar log for the abundance line plot.

    Returns ``None`` for non-positive / non-finite input so Plotly renders a gap
    (log of zero/negative is undefined), matching how the line plot already
    treats missing positions.
    """
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(v) or v <= 0:
        return None
    return float(_log_base_fn(log_base)(v))


def round_sig(x, sig=2):
    if x == 0:
        return 0
    exponent = int(np.floor(np.log10(abs(x))))
    factor = 10**(exponent - sig + 1)
    return round(x / factor) * factor


def round_sig_ceil(x, sig=2):
    if x == 0:
        return 0
    exponent = int(np.floor(np.log10(abs(x))))
    factor = 10**(exponent - sig + 1)
    return np.ceil(x / factor) * factor


def round_sig_floor(x, sig=2):
    if x == 0:
        return 0
    exponent = int(np.floor(np.log10(abs(x))))
    factor = 10**(exponent - sig + 1)
    return np.floor(x / factor) * factor


# Function to calculate y-ticks based on the min and max values of datasets
def calculate_y_ticks(min_values, max_values, log_transform, log_base=10):
    """
    Calculate y-ticks with proper handling of zero/negative values.

    ``log_base`` (10 or 2) selects the base for the log-scale tick anchors so
    the axis ladders on powers of 2 when Log2 is chosen, powers of 10 otherwise.
    """
    base = 2 if int(log_base) == 2 else 10
    # Check for empty lists or all zeros
    if not min_values or not max_values:
        return [0, 1, 10]  # Default scale if no data

    # Filter out zeros and negative values, keep valid positives
    valid_mins = [x for x in min_values if x > 0]
    valid_maxs = [x for x in max_values if x > 0]

    # If no valid values after filtering
    if not valid_mins or not valid_maxs:
        return [0, 1, 10]

    try:
        overall_min = np.nanmin(valid_mins)
        overall_max = np.nanmax(valid_maxs)

        # If min or max are zero or negative after nanmin/nanmax


        if log_transform:
            if overall_min <= 0 or overall_max <= 0:
                return [0, 1, 10]
            # Apply log transformation for log scale (base 10 or 2)
            log_b = np.log10 if base == 10 else np.log2
            min_power = base ** np.floor(log_b(overall_min))
            max_power = base ** np.ceil(log_b(overall_max))

            # Calculate midpoint in log space
            mid_point = np.sqrt(min_power * max_power)
            mid_point_rounded = base ** np.round(log_b(mid_point))

            return [min_power, mid_point_rounded, max_power]
        else:
            # For linear scale, return min, mid, max rounded to 2 sigfigs
            mid_point = (overall_min + overall_max) / 2
            overall_min_rounded = round_sig_floor(overall_min, 2)
            overall_max_rounded = round_sig(overall_max, 2)
            mid_point_rounded = round_sig_floor(mid_point, 2)
            return [overall_min_rounded, mid_point_rounded, overall_max_rounded]

    except (ValueError, RuntimeWarning):
        return [0, 1, 10]  # Fallback for any calculation errors


# Descending "nice" tick magnitudes within one decade, so a rounded axis bound
# reads cleanly (…, 6, 5, 4, 3, 2, 1 …). Used by calculate_signed_y_ticks.
_NICE_STEPS = (1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)


def calculate_signed_y_ticks(min_value, max_value):
    """Symmetric-about-zero y-ticks for the signed differential track.

    Unlike :func:`calculate_y_ticks` (positive-only abundance), SMD / log2FC
    values are small linear numbers straddling 0. Returns ``[-m, 0, m]`` where
    ``m`` is the max magnitude padded ~10% and rounded up to a clean tick, so the
    **zero baseline sits at the axis centre** and positive/negative excursions
    are directly comparable. Falls back to ``[-1, 0, 1]`` for empty/degenerate
    input.
    """
    vals = [v for v in (min_value, max_value)
            if v is not None and np.isfinite(v)]
    magnitude = max((abs(v) for v in vals), default=0.0)
    if magnitude == 0:
        return [-1.0, 0.0, 1.0]

    magnitude *= 1.1  # headroom so the extreme point isn't on the frame
    exp = np.floor(np.log10(magnitude))
    base = 10.0 ** exp
    mantissa = magnitude / base
    nice = next((s for s in _NICE_STEPS if mantissa <= s), 10.0)
    m = round(nice * base, 10)
    return [-m, 0.0, m]


_NON_AA_RE = re.compile(r'[^A-Za-z]')


def _resolve_index_range(start_val, end_val, sequence_length, position_offset=0):
    """
    Map a dataset peptide interval onto 0-based indices of the sequence as
    plotted, applying ``position_offset``.

    ``position_offset`` is the number of residues that were removed from the
    FRONT of the sequence *and* that the dataset's own numbering still
    counts. It is 0 when the dataset is already numbered from the mature
    protein (its position 1 == plotted residue 1) and ``strip_len`` when the
    dataset is numbered from the full precursor (its position 1 == the first
    residue of the signal peptide, which is no longer plotted). See
    :func:`detect_dataset_numbering`.

    Returns ``(start_idx, stop_idx)``, or ``(None, None)`` when the peptide
    lies entirely outside the plotted sequence. A peptide straddling the
    cleavage site is clipped to the part that survives.
    """
    start_idx = start_val - 1 - position_offset
    stop_idx = end_val - position_offset
    if stop_idx <= 0 or start_idx >= sequence_length or start_idx >= stop_idx:
        return None, None
    start_idx = max(start_idx, 0)
    stop_idx = min(stop_idx, sequence_length)
    if start_idx >= stop_idx:
        return None, None
    return start_idx, stop_idx


def detect_dataset_numbering(full_sequence, peptide_dataframe, strip_len):
    """
    Decide whether the dataset's peptide positions are numbered against the
    FULL precursor sequence or against the mature (signal-peptide-removed)
    sequence, so "Strip start sequence" can serve both real-world cases:

    * **'precursor'** — positions already line up with the untrimmed FASTA
      (the plot is correctly aligned *before* stripping). The user only wants
      the signal peptide gone to compact the figure, so every peptide must be
      renumbered down by ``strip_len`` to stay put → ``position_offset = strip_len``.
    * **'mature'** — positions are numbered from the mature protein while the
      FASTA still carries the signal peptide (the plot is misaligned *before*
      stripping). Trimming the sequence alone fixes it; positions must NOT
      move → ``position_offset = 0``.

    Detection is exact when the dataset carries a ``Sequence`` column: each
    peptide's own letters are compared against both hypotheses and the
    better-matching one wins. Without that column it falls back to bounds
    heuristics.

    Returns ``(numbering, reason)`` where numbering is 'precursor' or 'mature'.
    """
    trimmed = full_sequence[strip_len:]
    if strip_len <= 0:
        return 'mature', 'no residues stripped'
    if 'start' not in peptide_dataframe.columns or 'end' not in peptide_dataframe.columns:
        return 'mature', 'no start/end columns'

    # --- exact test: compare each peptide's own letters under both hypotheses
    if 'Sequence' in peptide_dataframe.columns:
        precursor_hits = mature_hits = 0
        for _, row in peptide_dataframe.iterrows():
            pep = row.get('Sequence')
            if not isinstance(pep, str) or not pep.strip():
                continue
            try:
                s = int(row['start'])
                e = int(row['end'])
            except (TypeError, ValueError):
                continue
            expected = _NON_AA_RE.sub('', pep).upper()
            if not expected:
                continue
            if 0 <= s - 1 and e <= len(full_sequence) \
                    and full_sequence[s - 1:e].upper() == expected:
                precursor_hits += 1
            if 0 <= s - 1 and e <= len(trimmed) \
                    and trimmed[s - 1:e].upper() == expected:
                mature_hits += 1
        if precursor_hits or mature_hits:
            if precursor_hits > mature_hits:
                return 'precursor', (
                    f'{precursor_hits} peptide sequence(s) match the full precursor numbering '
                    f'(vs {mature_hits} for mature numbering)')
            return 'mature', (
                f'{mature_hits} peptide sequence(s) match the mature-protein numbering '
                f'(vs {precursor_hits} for precursor numbering)')

    # --- fallback: no usable Sequence column, reason from the position bounds
    starts = pd.to_numeric(peptide_dataframe['start'], errors='coerce').dropna()
    ends = pd.to_numeric(peptide_dataframe['end'], errors='coerce').dropna()
    if ends.empty or starts.empty:
        return 'mature', 'no usable positions'
    if int(ends.max()) > len(trimmed):
        # Positions run past the end of the trimmed sequence, so they cannot
        # be numbered against it — they must be precursor coordinates.
        return 'precursor', (
            f'peptide positions run to {int(ends.max())}, past the {len(trimmed)}-residue '
            'mature sequence, so they must be numbered against the full precursor')
    if int(starts.min()) <= strip_len:
        # A peptide starting inside the (cleaved) signal-peptide region only
        # makes sense if positions are mature-numbered.
        return 'mature', (
            f'a peptide starts at position {int(starts.min())}, inside the {strip_len}-residue '
            'stripped region, so positions must be mature-numbered')
    return 'mature', 'ambiguous without a Sequence column — assumed mature-protein numbering'


def calculate_abundance(protein_sequence, peptide_dataframe, grouping_variable, position_offset=0):
    protein_sequence_length = len(protein_sequence)
    data = []
    col_names = []
    Avg_column = f'Avg_{grouping_variable}'

    for idx, row in peptide_dataframe.iterrows():
        try:
            start_val = int(row['start'])
            end_val = int(row['end'])
            # position_offset renumbers dataset positions onto protein_sequence
            # (see _resolve_index_range / detect_dataset_numbering).
            start_idx, stop_idx = _resolve_index_range(
                start_val, end_val, protein_sequence_length, position_offset)
            if start_idx is None:
                continue
            abundance_value = row[Avg_column]

            values = [0] * protein_sequence_length
            if abundance_value > 0:
                values[start_idx:stop_idx] = [abundance_value] * (stop_idx - start_idx)
            data.append(values)

            # Build column name in the same loop to guarantee alignment
            uid_suffix = ''
            if 'Unique Peptide ID' in row and not pd.isnull(row['Unique Peptide ID']):
                uid = str(row['Unique Peptide ID']).strip()
                if uid:
                    uid_suffix = f' {uid}'
            col_names.append(f"{start_val}-{end_val}{uid_suffix}")

        except (KeyError, ValueError) as e:
            error_message = f'{type(e).__name__}: {e} - Skipping this row'
            _log_validation(error_message)
            continue

    if not data:
        _log_validation("No valid data to process. Returning empty DataFrame.")
        return pd.DataFrame()

    # Create the initial DataFrame
    abundance_df = pd.DataFrame(data).T
    abundance_df.columns = col_names


    # Initialize count list
    count_list = []

    # Check if abundance_df is empty or has no non-zero values
    if abundance_df.empty or abundance_df.eq(0).all().all():
        # Create all zeros for counts and averages
        count_list = [0] * len(protein_sequence)
        abundance_df['average'] = 0
    else:
        # Calculate counts and averages normally
        count_list = abundance_df.gt(0).sum(axis=1)
        abundance_df['average'] = abundance_df.replace(0, np.nan).mean(axis=1)

    # Assign the counts and amino acids
    abundance_df['count'] = count_list
    abundance_df['AA'] = list(protein_sequence)
    return abundance_df


def _differs_only_by_il(expected, actual):
    """True when two equal-length residue strings differ *only* at positions
    where one side is I and the other L.

    Leucine and isoleucine are isobaric and cannot be told apart by mass
    spectrometry, so a search engine's I/L call is essentially a coin flip.
    An "expected vs actual" difference that is nothing but I<->L swaps is not
    a real misalignment and must not be flagged.
    """
    if len(expected) != len(actual):
        return False
    diffs = [(e, a) for e, a in zip(expected, actual) if e != a]
    return bool(diffs) and all({e, a} == {'I', 'L'} for e, a in diffs)


def validate_peptide_alignment(protein_sequence, peptide_dataframe, position_offset=0,
                               protein_id=None, max_examples=5):
    """
    Sanity-check that peptide start/end positions actually land inside
    ``protein_sequence`` (the sequence as it will be plotted — already
    trimmed by "Strip start sequence" if that setting is on), and, when the
    dataset carries the peptide's own ``Sequence`` text, that the residues
    at that position match it.

    ``position_offset`` is applied exactly as in the drawing functions (see
    :func:`_resolve_index_range`), so this validates what will actually be
    plotted rather than the raw dataset coordinates.

    This is the single place that flags a misaligned dataset (wrong manual
    offset, wrong UniProt signal-peptide length, or a FASTA that doesn't
    match the dataset's numbering at all) so users get one clear message
    instead of peptides silently disappearing from or misplaced on the plot.

    Returns a list with zero or one human-readable warning string (a single
    collapsible summary, not one line per peptide, so a systematic mismatch
    doesn't flood the UI). Differences that are purely I<->L swaps are
    ignored (see :func:`_differs_only_by_il`).
    """
    if 'start' not in peptide_dataframe.columns or 'end' not in peptide_dataframe.columns:
        return []

    protein_sequence_length = len(protein_sequence)
    has_seq_col = 'Sequence' in peptide_dataframe.columns

    out_of_range = 0
    mismatched = 0
    examples = []

    for _, row in peptide_dataframe.iterrows():
        try:
            start_val = int(row['start'])
            end_val = int(row['end'])
        except (TypeError, ValueError):
            continue

        raw_start_idx = start_val - 1 - position_offset
        start_idx, stop_idx = _resolve_index_range(
            start_val, end_val, protein_sequence_length, position_offset)

        if start_idx is None:
            if position_offset and end_val - position_offset <= 0:
                # Peptide lies entirely within the stripped signal peptide —
                # expected to disappear, not a misalignment.
                continue
            out_of_range += 1
            if len(examples) < max_examples:
                examples.append(
                    f"{start_val}-{end_val} falls outside the {protein_sequence_length}-"
                    "residue sequence"
                )
            continue

        if raw_start_idx < 0:
            # Straddles the cleavage site: only part of it survives the trim,
            # so its letters can't match the full peptide. Expected, skip.
            continue

        if has_seq_col:
            peptide_seq = row.get('Sequence')
            if isinstance(peptide_seq, str) and peptide_seq.strip():
                expected = _NON_AA_RE.sub('', peptide_seq).upper()
                actual = protein_sequence[start_idx:stop_idx].upper()
                if expected and expected != actual and not _differs_only_by_il(expected, actual):
                    mismatched += 1
                    if len(examples) < max_examples:
                        examples.append(f"{start_val}-{end_val} is '{actual}', expected '{expected}'")

    total = out_of_range + mismatched
    if not total:
        return []

    who = protein_id or "This protein"
    bullets = "".join(f"<li>{html.escape(e)}</li>" for e in examples)
    if total > len(examples):
        bullets += f"<li>+{total - len(examples)} more</li>"
    return [
        '<details class="hm-msg-details">'
        f"<summary>{html.escape(who)}: the dataset appears misaligned with the "
        "FASTA/UniProt sequence (click to expand)</summary>"
        f"<p>{total} peptide position(s) don&rsquo;t match the protein sequence:</p>"
        f"<ul>{bullets}</ul>"
        "<p>Is the protein sequence correct after the start/stop sequence has been "
        "stripped? Things to try:</p>"
        "<ul>"
        "<li>adjust the <strong>Strip Start Sequence</strong> value (or the UniProt "
        "signal-peptide length)</li>"
        "<li>upload a <strong>FASTA file</strong> whose sequence matches the "
        "numbering used in your dataset</li>"
        "</ul>"
        "</details>"
    ]


# ---------------------------------------------------------------------------
# Differential-comparison track (SMD / log2FC). See
# manuscript/docs/HEATMAP_DIFFERENTIAL_STATS.md.
# ---------------------------------------------------------------------------

def calculate_replicate_position_means(protein_sequence, peptide_dataframe, replicate_columns,
                                       position_offset=0):
    """Per-position mean abundance for each replicate column.

    For every replicate column, each peptide paints its abundance across the
    positions it spans; positions are then averaged over the covering peptides.
    A value of 0 means "not detected" and is treated as missing (``NaN``) so it
    is dropped from the mean — matching the single-series ``calculate_abundance``
    convention and the LOCKED "drop" missing-value policy.

    Returns a DataFrame indexed 0..len-1 (sequence position) with one column per
    replicate. Positions covered by no peptide (or only zero-abundance peptides)
    are ``NaN`` in that replicate and are dropped downstream by ``cohens_d`` /
    ``log2_fold_change``, so they never contribute a spurious zero.
    """
    seq_len = len(protein_sequence)
    cols = [c for c in (replicate_columns or []) if c in peptide_dataframe.columns]
    means = pd.DataFrame(index=range(seq_len), columns=cols, dtype=float)
    if not cols:
        return means

    for rep in cols:
        # position x peptide matrix for this replicate; 0 -> NaN so absent
        # peptides do not drag the per-position mean down.
        painted = []
        for _, row in peptide_dataframe.iterrows():
            try:
                start_idx, stop_idx = _resolve_index_range(
                    int(row['start']), int(row['end']), seq_len, position_offset)
            except (KeyError, ValueError, TypeError):
                continue
            if start_idx is None:
                # Falls outside the plotted sequence — already surfaced by
                # validate_peptide_alignment().
                continue
            val = pd.to_numeric(row.get(rep), errors='coerce')
            if not (pd.notna(val) and val > 0):
                continue
            column = np.full(seq_len, np.nan)
            column[start_idx:stop_idx] = val
            painted.append(column)
        if painted:
            # Uncovered positions are all-NaN columns; nanmean warns on those,
            # which is the intended "gap" outcome — silence the expected warning.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', category=RuntimeWarning)
                means[rep] = np.nanmean(np.vstack(painted), axis=0)
    return means


def calculate_contrast_track(means_a, means_b, metric='smd', eps=0.0):
    """Per-position signed contrast between two groups.

    ``means_a`` / ``means_b`` are the per-replicate position-mean frames from
    :func:`calculate_replicate_position_means` (they may have different replicate
    counts). Rows are aligned by integer sequence position; when the two frames
    differ in length (e.g. protein variants of unequal length) positions beyond
    the shorter frame are ``NaN``.

    ``metric='smd'`` → Cohen's *d* (the effect size); ``'log2fc'`` → log2 fold
    change. **Positive => higher in A.** Returns a 1-D ``numpy`` array of length
    ``max(len(a), len(b))``; undefined positions are ``NaN`` (rendered as gaps).
    """
    from data_analysis.services import stats as _stats

    n = max(len(means_a), len(means_b))
    out = np.full(n, np.nan)
    for p in range(n):
        a = means_a.iloc[p].to_numpy(dtype=float) if p < len(means_a) else np.array([])
        b = means_b.iloc[p].to_numpy(dtype=float) if p < len(means_b) else np.array([])
        if metric == 'log2fc':
            out[p] = _stats.log2_fold_change(a, b, eps=eps)
        else:
            out[p] = _stats.cohens_d(a, b)
    return out


def calculate_function(protein_sequence, peptide_dataframe, grouping_variable, position_offset=0):
    protein_sequence_length = len(protein_sequence)
    data = []
    col_names = []

    for _, row in peptide_dataframe.iterrows():
        start_val = int(row['start'])
        end_val = int(row['end'])
        start_idx, stop_idx = _resolve_index_range(
            start_val, end_val, protein_sequence_length, position_offset)
        if start_idx is None:
            # Falls outside the plotted sequence — already surfaced by
            # validate_peptide_alignment().
            continue

        function_value = row['function'] if 'function' in peptide_dataframe.columns else np.nan

        values = [None] * protein_sequence_length
        for i in range(start_idx, stop_idx):
            values[i] = function_value

        data.append(values)

        uid_suffix = ''
        if 'Unique Peptide ID' in row and not pd.isnull(row['Unique Peptide ID']):
            uid = str(row['Unique Peptide ID']).strip()
            if uid:
                uid_suffix = f' {uid}'
        col_names.append(f"{start_val}-{end_val}{uid_suffix}")

    function_df = pd.DataFrame(data).T
    function_df.columns = col_names
    return function_df

def export_heatmap_data_to_dict(protein_id, group_key, group_info, protein_sequence,
                               protein_species, protein_name, protein_df, is_all_null,
                               position_offset=0):
    """
    Exports the heatmap data to a dictionary based on filter type.

    protein_sequence is exactly what gets plotted — the full FASTA/UniProt
    sequence, or the "Strip start sequence" reduced sequence (front N
    residues physically removed by the caller) when that setting is on.
    position_offset renumbers protein_df's start/end onto that sequence: 0
    when the dataset is already numbered from the mature protein, strip_len
    when it is numbered from the full precursor. See detect_dataset_numbering.
    """
    grouping_var = group_info['grouping_variable']
    relevant_columns = group_info['abundance_columns']

    # Check if required columns exist, if not create default DataFrame
    required_columns = ['function', 'start', 'end']
    if not all(col in protein_df.columns for col in required_columns):
        filtered_df = pd.DataFrame({
            'start': [0],
            'end': [1],
            'function': ['0']
        })
        filtered_protein_df = filtered_df.copy()
    else:
        filtered_protein_df = protein_df.dropna(subset=['function', 'start', 'end'])
        if 'Unique Peptide ID' in filtered_protein_df.columns:
            filtered_df = filtered_protein_df[['start', 'end', 'function', 'Unique Peptide ID']]
        else:
            filtered_df = filtered_protein_df[['start', 'end', 'function']]

    # Calculate the heatmap data
    func_heatmap_df = calculate_function(
        protein_sequence, filtered_df, grouping_var, position_offset=position_offset)
    heatmap_df = calculate_abundance(
        protein_sequence, protein_df, grouping_var, position_offset=position_offset)
    filtered_heatmap_df = calculate_abundance(
        protein_sequence, filtered_protein_df, grouping_var, position_offset=position_offset)

    # Create the dictionary to store the results
    heatmap_data = {
        'protein_id': protein_id,
        'protein_sequence': protein_sequence,
        'protein_name': protein_name,
        'protein_species': protein_species,
        'func_heatmap_df': func_heatmap_df,
        'heatmap_df': heatmap_df,
        'filtered_heatmap_df': filtered_heatmap_df,
        # Raw per-peptide rows (start/end + any replicate 'Grouped:' columns) for
        # this protein. The differential-comparison track needs these to build
        # per-replicate per-position means; the position-painted heatmap_df above
        # only carries the Avg_ column, not the individual replicates.
        'peptide_df': protein_df,
    }

    return heatmap_data

def process_available_data(available_data_variables_dict, filter_type, selected_functions, selected_peptides, bio_or_pep, log_transform, manual_y_axis=False, y_min_manual=None, y_max_manual=None, log_base=10):
    """
    Process data for visualization when update_plot is called
    Returns: dict with data and 'errors' key containing list of error messages
    """
    errors = []
    if not available_data_variables_dict:
        return {'errors': ["No data available for processing"]}

    count_list = []
    min_values = []
    max_values = []
    seq_len_list = []

    # Process each variable's data for visualization
    for var in available_data_variables_dict:
        # Default to the unfiltered heatmap so an unrecognised filter_type falls back
        # to showing all peptides rather than raising UnboundLocalError below.
        df = available_data_variables_dict[var]['heatmap_df']
        if filter_type == 'all-peptides':
            df = available_data_variables_dict[var]['heatmap_df']
        elif filter_type == 'bioactive-only':
            df = available_data_variables_dict[var]['filtered_heatmap_df']
        elif filter_type == 'functional-only':
            # Filter both functional and Abundance data for specific functions
            func_df = available_data_variables_dict[var]['function_heatmap_df']
            if func_df is not None and not func_df.empty:
                # Get list of columns that contain any of the selected functions
                functional_positions = []
                for col in func_df.columns:
                    # Check if any selected function appears in any cell of this column
                    if any(any(func in str(cell) for func in selected_functions) for cell in func_df[col]):
                        functional_positions.append(col)

                # Get Abundance data for these positions
                abs_df = available_data_variables_dict[var]['filtered_heatmap_df']
                if abs_df is not None and not abs_df.empty and functional_positions:
                    # Only keep columns that exist in abs_df
                    valid_columns = [col for col in functional_positions if col in abs_df.columns]
                    if valid_columns:
                        # Get the filtered Abundance data for selected columns
                        selected_abs_df = abs_df[valid_columns]

                        # Create DataFrame with sequence and selected columns
                        df = pd.DataFrame({
                            'AA': abs_df['AA'],
                        })

                        # Add the selected columns
                        df = pd.concat([df, selected_abs_df], axis=1)

                        # Now calculate count and average from all data
                        non_zero_mask = df.drop('AA', axis=1) > 0
                        df['count'] = non_zero_mask.sum(axis=1)
                        df['average'] = df.drop('AA', axis=1).where(non_zero_mask).mean(axis=1)
                    else:
                        df = None
                else:
                    df = None
            else:
                df = None

        # When specific peptide intervals are selected, narrow the heatmap data to only
        # those intervals so that peptide_counts and ms_data (heatmap colours + line
        # averages) reflect the selection rather than all peptides.
        if bio_or_pep == '1' and selected_peptides and df is not None and not df.empty:
            try:
                pep_cols = [col for col in selected_peptides if col in df.columns]
                if pep_cols:
                    keep = [c for c in ['AA'] + pep_cols if c in df.columns]
                    df = df[keep].copy()
                    non_zero = df[pep_cols].gt(0)
                    df['count']   = non_zero.sum(axis=1)
                    df['average'] = df[pep_cols].where(non_zero).mean(axis=1)
            except Exception:
                pass  # fall back to the unfiltered dataframe on any error

        # Get counts and MS data
        if df is None or df.empty or 'count' not in df.columns or 'average' not in df.columns:
            # Create default DataFrame with zero values
            protein_length = len(available_data_variables_dict[var]['protein_sequence'])
            df = pd.DataFrame({
                'count': [0] * protein_length,
                'average': [0] * protein_length,
                'AA': list(available_data_variables_dict[var]['protein_sequence'])
            })

        # Get counts and MS data
        peptide_counts = df['count']
        ms_data = df['average']
        count_list.append(peptide_counts)

        # Calculate min/max MS values
        min_ms = ms_data[ms_data > 0].min()
        max_ms = ms_data.max()
        min_values.append(min_ms)
        max_values.append(max_ms)

        # Add computed properties to the variable data
        available_data_variables_dict[var].update({
            'peptide_counts': peptide_counts,
            'ms_data': ms_data,
            'ms_data_list': list(ms_data),
            'AA_list': df['AA'].tolist(),
            'max_peptide_counts': peptide_counts.max(),
            'min_peptide_counts': peptide_counts.min(),
            'max_ms_data': max_ms,
            'min_ms_data': min_ms,
        })

        # Process bioactive peptide data
        columns_to_include = df.columns.difference(['AA', 'count'])
        df_filtered = df[columns_to_include]

        available_data_variables_dict[var].update({
            'bioactive_peptide_abs_df': df_filtered,
            'bioactive_peptide_func_df': available_data_variables_dict[var]['function_heatmap_df']
        })

        seq_len_list.append(len(available_data_variables_dict[var]['AA_list']))

    # Calculate global values
    global y_ticks
    max_sequence_length = max(seq_len_list)

    # Process counts
    if len(set([item for sublist in count_list for item in sublist])) >= 1:
        flat_list = [item for sublist in count_list for item in sublist]
        # Check if there are any non-zero values before filtering
        non_zero_values = [item for item in flat_list if item != 0]
        if non_zero_values:  # If there are non-zero values
            list_of_counts = set(non_zero_values)
            max_count = max(non_zero_values)
            num_unique_count = len(set(non_zero_values))
        else:  # If all values are zero
            list_of_counts = {0}
            max_count = 0
            num_unique_count = 1
    else:
        max_count = 0
        num_unique_count = 0
        list_of_counts = set()
    # Calculate number of colors
    num_colors = 6 if num_unique_count >= 6 else num_unique_count

    # Override min/max values with manual settings if enabled
    if manual_y_axis:
        min_values = [y_min_manual]
        max_values = [y_max_manual]
    # Calculate global values
    if min_values and max_values:
        y_ticks = calculate_y_ticks(min_values, max_values, log_transform, log_base)
        y_ticks_str = ', '.join(f'{tick:.2e}' for tick in y_ticks)
        y_ticks_html = f'<b>Max/Min of MS data (y-ticks):</b> {y_ticks_str}'
    else:
        y_ticks = [0.1, 1, 10]  # Default log scale values
        y_ticks_html = '<span style="color:red;">Insufficient data to calculate MS data y-ticks.</span>'
        y_ticks_html += f'<b>Protein Sequence Length:</b> {max_sequence_length}'


    return {
        'list_of_counts': list_of_counts,
        'min_values': min_values,
        'max_values': max_values,
        'seq_len_list': seq_len_list,
        'max_sequence_length': max_sequence_length,
        'y_ticks': y_ticks,
        'max_count': max_count,
        'num_unique_count': num_unique_count,
        'num_colors': num_colors,
        'errors': errors

    }

# ---------------------------------------------------------------------------
# Helper: derive a short protein label from the variable dict.
# Mirrors the notebook's `process_data_variables` / `protein_name_short` logic.
# Used for the x-axis label ("<name> Sequence"). Protein names are already
# cleaned at load time (data_processor._clean_protein_name populates
# protein_dict with the short name), so nothing is stripped here.
# ---------------------------------------------------------------------------
def _derive_protein_short_name(available_data_variables_dict: dict) -> str:
    """Return the short protein label (e.g. 'Beta-casein') from the variable
    dict's protein_name / protein_id fields, or '' if none.

    Rules (matching the notebook): multiple IDs sharing one name → that name;
    multiple IDs with distinct names → names joined by '_'; otherwise the sole
    name, else the sole protein ID.
    """
    pid_list   = []
    pname_list = []
    for vd in available_data_variables_dict.values():
        pid   = vd.get('protein_id',   '')
        pname = str(vd.get('protein_name', '')).strip()
        if pid:
            pid_list.append(str(pid))
        if pname:
            pname_list.append(pname)

    unique_ids   = list(dict.fromkeys(pid_list))
    unique_names = list(dict.fromkeys(pname_list))

    if len(unique_ids) > 1 and len(unique_names) == 1:
        return unique_names[0]
    if len(unique_ids) > 1 and len(unique_names) > 1:
        return '_'.join(unique_names)
    if unique_names:
        return unique_names[0]
    if unique_ids:
        return unique_ids[0]
    return ''


def _derive_xaxis_label(available_data_variables_dict: dict) -> str:
    """Return '<ProteinName> Sequence' derived from protein_name / protein_id fields."""
    short = _derive_protein_short_name(available_data_variables_dict)
    return f"{short} Sequence" if short else 'Protein Sequence Position'


# Update plot function
def update_plot(available_data_variables_dict, ms_average_choice, bio_or_pep, selected_peptides, selected_functions,
                hm_selected_color, lp_selected_color, avglp_selected_color,
                xaxis_label, yaxis_label, yaxis_position, legend_title_input_1, legend_title_input_2,
                legend_title_input_3, filter_type, log_transform,
                manual_y_axis, y_min_manual, y_max_manual, log_base=10, plot_compact=False,
                plot_landscape_interactive=False, plot_portrait_interactive=False,
                plot_title='', aa_on_tiles=False,
                comparison_mode=False, series_a=None, series_b=None,
                comparison_metric='smd',
                font_size_xaxis_label=14, font_size_yaxis_label=14,
                font_size_legend=14, font_size_plot_title=16,
                font_size_xaxis_tick=14, font_size_yaxis_tick=14,
                font_size_var_label=14,
                legend_position='right'):
    # Every heatmap orientation is now rendered by Plotly
    # (visualize_sequence_heatmap_interactive / _compact); the matplotlib
    # landscape/portrait renderers and their chunking machinery were retired.
    import matplotlib.pyplot as plt
    all_errors = []

    if not available_data_variables_dict:  # Check if data is available
        all_errors.append('No data is available for plotting.')
        return (None, None, None, None, None, all_errors, [])  # Return tuple with errors

    # If the caller supplied no x-axis label, auto-derive it from protein name(s)
    # exactly as the notebook does when it sets xaxis_label_input to
    # f"{protein_name_short} Sequence".
    if not (xaxis_label or '').strip():
        xaxis_label = _derive_xaxis_label(available_data_variables_dict)

    result = process_available_data(available_data_variables_dict, filter_type, selected_functions,  selected_peptides, bio_or_pep, log_transform, manual_y_axis, y_min_manual, y_max_manual, log_base=log_base)
    # Interactive Plotly returns (compact / portrait / landscape). The two
    # leading None slots preserve the historical 7-tuple shape (they used to
    # carry the matplotlib portrait/landscape figures, now always None).
    fig_compact_return      = None
    fig_port_plotly_return  = None
    fig_land_plotly_return  = None

    # Unpack the dictionary into individual global variables
    if result:
        # Check for processing errors
        if 'errors' in result and result['errors']:
            all_errors.extend(result['errors'])

        # Unpack the dictionary into individual global variables
        global list_of_counts, min_values, max_values, max_sequence_length
        global max_count, num_unique_count, num_colors, y_ticks

        list_of_counts = result['list_of_counts']
        min_values = result['min_values']
        max_values = result['max_values']
        max_sequence_length = result['max_sequence_length']
        y_ticks = result['y_ticks']
        max_count = result['max_count']
        num_unique_count = result['num_unique_count']
        num_colors = result['num_colors']

        # ── Differential-comparison track (SMD / log2FC) ─────────────────────
        # When comparison mode is on and both series are replicate-backed, build
        # a single signed contrast line (per-position) that REPLACES the abundance
        # line in the row-1 panel. Signed data → symmetric-about-0 axis, zero
        # baseline, and NO log-transform (guarded below). Any validation failure
        # leaves comparison off and the normal abundance path untouched.
        comparison_track = None
        comparison_render_dict = None
        if comparison_mode:
            metric = comparison_metric if comparison_metric in ('smd', 'log2fc') else 'smd'
            va = available_data_variables_dict.get(series_a)
            vb = available_data_variables_dict.get(series_b)
            if not series_a or not series_b or series_a == series_b:
                all_errors.append('Comparison requires two distinct series (Series A and Series B).')
            elif va is None or vb is None:
                all_errors.append('Comparison series not found in the current protein/variable selection.')
            elif not va.get('replicate_columns') or not vb.get('replicate_columns'):
                all_errors.append('Comparison needs replicate-level (Grouped) data for both series; this looks like a single-average file. Re-export from Data Transformation with grouped replicate columns.')
            elif va.get('peptide_df') is None or vb.get('peptide_df') is None:
                all_errors.append('Comparison unavailable: raw peptide data missing for one of the selected series.')
            else:
                try:
                    means_a = calculate_replicate_position_means(
                        va['protein_sequence'], va['peptide_df'], va['replicate_columns'],
                        position_offset=va.get('position_offset', 0))
                    means_b = calculate_replicate_position_means(
                        vb['protein_sequence'], vb['peptide_df'], vb['replicate_columns'],
                        position_offset=vb.get('position_offset', 0))
                    track = calculate_contrast_track(means_a, means_b, metric=metric)
                    finite = track[np.isfinite(track)]
                    if finite.size == 0:
                        all_errors.append('Comparison produced no defined positions (no overlapping replicate coverage).')
                    else:
                        lo = float(np.min(finite))
                        hi = float(np.max(finite))
                        y_ticks = calculate_signed_y_ticks(lo, hi)
                        log_transform = False  # signed data — never log-transform

                        # Series labels. When A and B are DIFFERENT proteins the
                        # bare variable name ("Bitter") does not distinguish the two
                        # density rows, so prefix the protein name; same protein →
                        # the variable name alone is unambiguous.
                        # Use the bare sample name (var_label) — not the display
                        # 'label', which is already protein-prefixed in multi-protein
                        # selections and would otherwise double up below.
                        var_a = va.get('var_label') or va.get('label', series_a)
                        var_b = vb.get('var_label') or vb.get('label', series_b)
                        prot_a = str(va.get('protein_name') or va.get('protein_id') or '').strip()
                        prot_b = str(vb.get('protein_name') or vb.get('protein_id') or '').strip()
                        if va.get('protein_id') != vb.get('protein_id'):
                            if var_a == var_b:
                                label_a, label_b = (prot_a or var_a), (prot_b or var_b)
                            else:
                                label_a = f'{prot_a} {var_a}'.strip()
                                label_b = f'{prot_b} {var_b}'.strip()
                        else:
                            label_a, label_b = var_a, var_b

                        comparison_track = {
                            'values': track,
                            'metric': metric,
                            'label_a': label_a,
                            'label_b': label_b,
                            'aa_list': list(va.get('AA_list') or vb.get('AA_list') or va['protein_sequence']),
                        }
                        # Restrict the density heatmap to ONLY the two compared
                        # series (not every selected variable), and carry the
                        # protein-aware labels onto the tile rows. Shallow-copy so
                        # the shared dict entries are not mutated.
                        a_entry = dict(va); a_entry['label'] = label_a
                        b_entry = dict(vb); b_entry['label'] = label_b
                        comparison_render_dict = {series_a: a_entry, series_b: b_entry}
                except Exception as exc:
                    all_errors.append(f'Error computing comparison track: {exc}')

        selected_legend_title = [legend_title_input_1, legend_title_input_2, legend_title_input_3, legend_title[4]]

        # Your plotting code here, using the widget values as inputs

        cmap = plt.get_cmap(hm_selected_color)
        avg_cmap = plt.get_cmap(avglp_selected_color)

        # Comprehensive validation of plotting configuration.
        # When a comparison track is active it supplies the entire row-1 line
        # panel (the abundance / individual-peptide lines are suppressed), so the
        # abundance-source scenarios below do not apply — skip them.
        config_errors = []

        # Scenario 1: No averaged data + no specific peptides selected.
        # filter_type='bioactive-only' is itself a valid data source (all bioactive
        # peptides), so it must be excluded from this check even when bio_or_pep='no'.
        if comparison_track is not None:
            pass
        elif ms_average_choice == 'no' and (not bio_or_pep or bio_or_pep == 'no') and filter_type not in ('bioactive-only', 'functional-only'):
            config_errors.append('Invalid configuration: "Plot Averaged Data" is set to "No" and "Plot Specific Peptides" is not selected. Please enable at least one data source for plotting.')

        # Scenario 2: Peptide intervals selected but no actual intervals provided
        elif bio_or_pep == '1' and (not selected_peptides or len(selected_peptides) == 0):
            config_errors.append('Invalid configuration: "Plot Specific Peptides" is set to "Peptide Intervals" but no peptide intervals are selected. Please select peptide intervals from the "Specific Options" dropdown or change the plotting configuration.')

        # Scenario 3: Bioactive functions selected but no actual functions provided
        elif bio_or_pep == '2' and (not selected_functions or len(selected_functions) == 0):
            config_errors.append('Invalid configuration: "Plot Specific Peptides" is set to "Bioactive Functions" but no bioactive functions are selected. Please select bioactive functions from the "Specific Options" dropdown or change the plotting configuration.')


        if config_errors:
            all_errors.extend(config_errors)
            return None, None, None, None, None, all_errors, []

        else:
            """ #                                       Function Call to Generate Plot
            #---------------------------------------------------------------------------------------------------------------------------------------------
            """

            # (matplotlib portrait/landscape rendering removed — Plotly only)

            # ── Interactive Plotly plots for portrait / landscape ─────────────
            _interactive_params = dict(
                xaxis_label      = xaxis_label,
                yaxis_label      = yaxis_label,
                legend_title     = selected_legend_title,
                cmap             = cmap,
                avg_cmap         = avg_cmap,
                lp_selected_color= lp_selected_color,
                selected_peptides= selected_peptides,
                selected_functions= selected_functions,
                ms_average_choice= ms_average_choice,
                bio_or_pep       = bio_or_pep,
                filter_type      = filter_type,
                log_transform    = log_transform,
                log_base         = log_base,
                y_ticks          = y_ticks,
                plot_title       = plot_title,
                aa_on_tiles      = aa_on_tiles,
                comparison_track = comparison_track,
                font_size_xaxis_label   = font_size_xaxis_label,
                font_size_yaxis_label   = font_size_yaxis_label,
                font_size_legend        = font_size_legend,
                font_size_plot_title    = font_size_plot_title,
                font_size_xaxis_tick    = font_size_xaxis_tick,
                font_size_yaxis_tick    = font_size_yaxis_tick,
                font_size_var_label     = font_size_var_label,
                legend_position         = legend_position,
            )

            # In comparison mode the density heatmap should show ONLY the two
            # compared series (with protein-aware labels), not every selected
            # variable — so render from the filtered dict when it was built.
            interactive_dict = comparison_render_dict if comparison_render_dict is not None else available_data_variables_dict

            fig_port_plotly = None
            if plot_portrait_interactive:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', UserWarning)
                    try:
                        fig_port_plotly = visualize_sequence_heatmap_interactive(
                            interactive_dict,
                            'Portrait',
                            **_interactive_params,
                        )
                    except Exception as exc:
                        all_errors.append(f'Error generating portrait interactive plot: {exc}')

            fig_land_plotly = None
            if plot_landscape_interactive:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', UserWarning)
                    try:
                        fig_land_plotly = visualize_sequence_heatmap_interactive(
                            interactive_dict,
                            'Landscape',
                            **_interactive_params,
                        )
                    except Exception as exc:
                        all_errors.append(f'Error generating landscape interactive plot: {exc}')

            # ── Compact interactive Plotly plot ───────────────────────────────
            fig_compact = None
            if plot_compact:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', UserWarning)
                    try:
                        fig_compact = visualize_sequence_heatmap_compact(
                            available_data_variables_dict,
                            xaxis_label,
                            yaxis_label,
                            selected_legend_title,
                            cmap,
                            y_ticks,
                            log_transform,
                            log_base,
                            plot_title=plot_title,
                            font_size_xaxis_label=font_size_xaxis_label,
                            font_size_yaxis_label=font_size_yaxis_label,
                            font_size_legend=font_size_legend,
                            font_size_plot_title=font_size_plot_title,
                            font_size_xaxis_tick=font_size_xaxis_tick,
                            legend_position=legend_position,
                        )
                    except Exception as exc:
                        all_errors.append(f'Error generating compact plot: {exc}')
                        fig_compact = None

        # Collect the interactive Plotly figures for return.
        if plot_portrait_interactive:
            fig_port_plotly_return = fig_port_plotly  # may be None if generation failed
        if plot_landscape_interactive:
            fig_land_plotly_return = fig_land_plotly  # may be None if generation failed
        if plot_compact:
            fig_compact_return = fig_compact

        if not (plot_compact or plot_landscape_interactive or plot_portrait_interactive):
            all_errors.append("Select at least one plot type to generate a plot.")

        # Collect missing data notifications for user feedback
        missing_data_notifications = collect_missing_data_notifications(
            available_data_variables_dict, selected_peptides, selected_functions, bio_or_pep
        )

        # 7-tuple shape kept for callers; first two slots (once matplotlib
        # portrait/landscape) are always None now.
        return (None, None, fig_compact_return,
                fig_port_plotly_return, fig_land_plotly_return,
                all_errors, missing_data_notifications)
    else:
        all_errors.append('Reselect "Variable Options", no data is available for plotting.')
        return None, None, None, None, None, all_errors, []  # Return None values with errors if no result

"""_________________________________________Data Visualization Functions_________________________________"""

# Function to extract non-zero, non-NaN values
def extract_non_zero_non_nan_values(df):
    unique_functions = set()
    # Iterate over each value in the DataFrame
    for value in df.stack().values:  # df.stack() stacks the DataFrame into a Series
        if value != 0 and not pd.isna(value):  # Check if value is non-zero and not NaN
            if isinstance(value, str):
                # If the value is a string, it could contain multiple delimited entries
                entries = value.split('; ')
                unique_functions.update(entries)
            else:
                unique_functions.add(value)
    return unique_functions

    # This function is used in plotting to filter the data to only plot selected peptides or bioactive peptides, independent of averaged data


def safe_str(value):
    """Safely convert any value to string, handling None, NaN, and Series objects"""
    try:
        if value is None:
            return ""

        # For Series objects, just return empty string to avoid truth value errors
        if hasattr(value, 'iloc'):
            return ""

        # Convert to string
        return str(value)
    except:
        return ""


def filter_data_by_selection(bp_abs, bp_func, selected_peptides, selected_functions, bio_or_pep, filter_type, ms_average_choice, var_name):
    """
    Final simplified version that uses a completely different approach to avoid Series truth value issues.
    Returns: (filtered_bp_abs, filtered_bp_func, errors)
    """
    errors = []

    # ---------- Handle peptide intervals selection (bio_or_pep = '1') ----------
    if bio_or_pep == '1' and selected_peptides:
        # Find our selected peptide columns that exist in the DataFrame
        if bp_abs is None or not isinstance(bp_abs, pd.DataFrame) or bp_abs.empty:
            errors.append(f"Abundance DataFrame is empty for variable: {var_name}")
            return pd.DataFrame(), pd.DataFrame(), errors

        peptide_columns = [col for col in selected_peptides if col in bp_abs.columns]
        non_matching_peptides = [col for col in selected_peptides if col not in bp_abs.columns]

        if peptide_columns:
            # Add metadata columns if they exist
            meta_columns = ['AA', 'count', 'average']
            keep_columns = [col for col in meta_columns if col in bp_abs.columns]
            all_columns = keep_columns + peptide_columns

            # Create new DataFrames with only the columns we want
            filtered_bp_abs = bp_abs[all_columns].copy()
            filtered_bp_func = pd.DataFrame()

            return filtered_bp_abs, filtered_bp_func, errors
        else:
            # No selected peptides exist in this protein/variable combination — silently skip.
            # This is expected when plotting multiple proteins and a peptide belongs to only one.
            return pd.DataFrame(), pd.DataFrame(), errors

    # ---------- Handle bioactive functions selection (bio_or_pep = '2') ----------
    elif bio_or_pep == '2' and selected_functions:
        if bp_func is None or not isinstance(bp_func, pd.DataFrame) or bp_func.empty:
            errors.append(f"Function DataFrame is empty for variable: {var_name}")
            return pd.DataFrame(), pd.DataFrame(), errors
        try:
            # Simple string-based approach that completely avoids pandas operations
            # This manually scans the DataFrame for columns containing the selected functions
            columns_with_functions = []
            functions_found = set()

            # Get list of column names first
            all_columns = list(bp_func.columns)

            # Iterate through the function DataFrame columns
            for col_name in all_columns:
                # Process a column at a time - get non-NA function values
                func_values = bp_func[col_name].dropna()

                if func_values.empty:
                    continue

                # Track matched functions for this column
                matched_functions = set()

                # Check for matches with selected functions, splitting on ';'
                for val in func_values:
                    # Convert to string (safely)
                    try:
                        value_str = str(val)
                    except:
                        continue

                    # Skip empty strings or 'nan'
                    if not value_str or value_str == 'nan':
                        continue

                    # Split the value on ';' and strip whitespace
                    sub_vals = [s.strip() for s in value_str.split(';')]

                    for sub_val in sub_vals:
                        # Only match if sub_val matches exactly one of the selected_functions (case-insensitive, ignoring whitespace)
                        for func in selected_functions:
                            if func.lower().replace(" ", "") == sub_val.lower().replace(" ", ""):
                                matched_functions.add(func)

                # If we found any matches in this column, add it to the list
                if matched_functions:
                    columns_with_functions.append(col_name)

                # Add matched functions to the global set
                functions_found.update(matched_functions)

            # Remove any duplicates
            columns_with_functions = list(set(columns_with_functions))

            # Identify non-matching functions
            non_matching_functions = [func for func in selected_functions if func not in functions_found]

            if columns_with_functions:
                # Only include columns that are in both DataFrames
                valid_cols_abs = [col for col in columns_with_functions if col in bp_abs.columns]
                valid_cols_func = [col for col in columns_with_functions if col in bp_func.columns]

                # Create filtered DataFrames
                if valid_cols_abs:
                    filtered_bp_abs = bp_abs[valid_cols_abs].copy()
                else:
                    filtered_bp_abs = pd.DataFrame()

                if valid_cols_func:
                    filtered_bp_func = bp_func[valid_cols_func].copy()
                else:
                    filtered_bp_func = pd.DataFrame()

                # Add warning for non-matching functions if any
                if non_matching_functions:
                    errors.append(f'Some selected functions not found in data for variable: {var_name}. Missing: {non_matching_functions}')

                return filtered_bp_abs, filtered_bp_func, errors
            else:
                errors.append(f'No columns found containing any of the selected functions for variable: {var_name}. Non-matching functions: {non_matching_functions}')
                return pd.DataFrame(), pd.DataFrame(), errors

        except Exception as e:
            errors.append(f"Error in function filtering for variable: {var_name}: {str(e)}")
            return pd.DataFrame(), pd.DataFrame(), errors

    else:
        return bp_abs, bp_func, errors

def get_grouped_colors(counts, max_count, num_groups, plot_zero, cmap):
    from matplotlib.colors import Normalize
    # Initialize colors list with None to maintain length
    colors = [None] * len(counts)

    # Set the start point based on the user's input
    start_point = 0 if plot_zero == 'yes' else 1

    # Group counts into fewer categories if necessary, excluding zeros
    group_bounds = np.linspace(start_point, max_count, num_groups + 1)
    group_labels = np.digitize(counts, bins=group_bounds, right=True)  # Find group labels for counts
    ##print("DEBUG: group_labels", group_labels)
    norm = Normalize(vmin=start_point, vmax=num_groups)  # Normalize based on the number of groups

    # Map each count to a color
    #for i, count in enumerate(counts):
    for i, (count, label) in enumerate(zip(counts, group_labels)):
        if plot_zero == 'no' and count == 0:
            colors[i] = 'white'
        else:
            # Find which group the count belongs to
            group_idx = np.digitize(count, group_bounds) - 1

            if max_count > 20:
                # For large ranges, use continuous color mapping
                colors[i] = cmap(count / max_count)
            else:
                # For smaller ranges, use discrete color groups
                #colors[i] = cmap(norm(group_bounds[min(group_idx, len(group_bounds)-1)]))
                colors[i] = cmap(norm(label))

    #print("DEBUG: label",label)
    return colors

# Function to create legend for heatmap
def create_heatmap_legend_handles(cmap, num_colors, max_count, plot_zero):
    """
    Create color-coded legend handles for the heatmap with proper handling of zero values
    and error cases
    """
    import matplotlib.patches as patches
    from matplotlib.colors import Normalize
    try:
        # Handle case where all values are zero
        if max_count == 0:
           norm = Normalize(vmin=0, vmax=1)
           color = cmap(norm(0))
           return [patches.Patch(color=color, label='0')], ['0']

        # Set the start point based on plot_zero
        start_point = 0 if plot_zero == 'yes' else 1

        # Create group boundaries to match get_grouped_colors
        count_ranges = np.linspace(start_point, max_count, num_colors + 1)

        # Determine interval type based on count vs intervals
        if max_count > num_colors:
           plt_interval = max_count
        else:
           plt_interval = num_colors

        legend_handles = []
        heatmap_labels = []
        norm = Normalize(vmin=start_point, vmax=max_count)

        # For small ranges, start from index 1 to skip the duplicate 1
        start_idx = 1 if plt_interval <= 6 else 0
        for i in range(start_idx, len(count_ranges)):
            color = cmap(norm(count_ranges[i]))
            if plt_interval <= 6 and plot_zero == 'no':
                label = f'{int(count_ranges[i])}'
            elif plt_interval <= 6 and plot_zero == 'yes':
                label = f'{int(count_ranges[i])}'
            elif i + 1 >= len(count_ranges):
                label = f'{int(count_ranges[i])} - {max(count_ranges)}'
                break
            else:
                # Create label showing range
                start_val = int(count_ranges[i])
                end_val = int(count_ranges[i + 1])
                label = f'{start_val} - {end_val}'

            #if i = len(count_ranges):
            #    color = cmap(norm(count_ranges[i]))

            #print("DEBUG: i:",i,"label:",label)
            #legend_handles.append(patches.Patch(color=color, label=label))
            legend_handles.append(patches.Patch(color=color, label=label))
            heatmap_labels.append(label)
        return legend_handles, heatmap_labels

    except Exception as e:
       # Fallback for error cases
       color = cmap(0)
       handle = patches.Patch(color=color, label='0')
       return [handle], ['0']


def collect_missing_data_notifications(available_data_variables_dict, selected_peptides, selected_functions, bio_or_pep):
    """
    Collect notifications about missing data for specific variables during plotting.
    Returns a list of notification messages for user display.
    """
    notifications = []

    if bio_or_pep == '1' and selected_peptides:  # Peptide intervals selected
        # Track which peptides are missing from which variables
        missing_by_variable = {}

        for var_name, var_data in available_data_variables_dict.items():
            missing_peptides = []
            if 'heatmap_df' in var_data:
                available_columns = var_data['heatmap_df'].columns.tolist()
                for peptide in selected_peptides:
                    if peptide not in available_columns:
                        missing_peptides.append(peptide)

            if missing_peptides:
                missing_by_variable[var_name] = missing_peptides

        # Create notification messages
        for var_name, missing_peptides in missing_by_variable.items():
            if len(missing_peptides) == 1:
                notifications.append(f"⚠️ Variable '{var_name}': Peptide interval '{missing_peptides[0]}' is not available in this dataset.")
            else:
                peptide_list = "', '".join(missing_peptides)
                notifications.append(f"⚠️ Variable '{var_name}': Peptide intervals '{peptide_list}' are not available in this dataset.")

    elif bio_or_pep == '2' and selected_functions:  # Bioactive functions selected
        # Track which functions are missing from which variables
        missing_by_variable = {}

        for var_name, var_data in available_data_variables_dict.items():
            missing_functions = []
            if 'function_heatmap_df' in var_data:
                func_df = var_data['function_heatmap_df']
                if func_df is not None and not func_df.empty:
                    # Extract all functions present in this variable
                    present_functions = set()
                    for col in func_df.columns:
                        for cell in func_df[col]:
                            if pd.notna(cell) and str(cell) != '0':
                                cell_str = str(cell)
                                if ';' in cell_str:
                                    present_functions.update(cell_str.split(';'))
                                else:
                                    present_functions.add(cell_str)

                    # Check which selected functions are missing
                    for function in selected_functions:
                        if function not in present_functions:
                            missing_functions.append(function)
            else:
                # If no function data at all, all selected functions are missing
                missing_functions = selected_functions.copy()

            if missing_functions:
                missing_by_variable[var_name] = missing_functions

        # Create notification messages
        for var_name, missing_functions in missing_by_variable.items():
            if len(missing_functions) == 1:
                notifications.append(f"⚠️ Variable '{var_name}': Bioactive function '{missing_functions[0]}' is not available in this dataset.")
            else:
                function_list = "', '".join(missing_functions)
                notifications.append(f"⚠️ Variable '{var_name}': Bioactive functions '{function_list}' are not available in this dataset.")

    return notifications


"""_________________________________________Landscape Plot________________________________________"""


"""_______________________________Interactive Plotly Plot (Landscape / Portrait)_____________"""
def visualize_sequence_heatmap_interactive(
    available_data_variables_dict,
    orientation,          # 'Landscape' or 'Portrait'
    xaxis_label,
    yaxis_label,
    legend_title,         # list [title1, title2, title3, ...]
    cmap,                 # matplotlib colormap for heatmap (peptide count colours)
    avg_cmap,             # matplotlib colormap for averaged line colours
    lp_selected_color,    # colourmap name for individual peptide lines
    selected_peptides,
    selected_functions,
    ms_average_choice,    # 'yes' | 'no' | 'only'
    bio_or_pep,           # 'no' | '1' | '2'
    filter_type,
    log_transform,
    y_ticks,
    log_base=10,          # base for the abundance log transform (10 or 2)
    plot_title='',        # optional user figure title; blank → no title
    aa_on_tiles=False,    # print the AA letter on each sample's colored tile
                          # (instead of the shared grey sequence strip) — lets you
                          # compare two proteins / variants row-by-row
    comparison_track=None,# dict(values, metric, label_a, label_b, aa_list) — when
                          # set, the row-1 panel shows a single signed SMD/log2FC
                          # contrast line (zero baseline, symmetric axis) INSTEAD
                          # of the abundance lines. See HEATMAP_DIFFERENTIAL_STATS.md.
    # Appearance Settings font-size overrides. font_size_legend sizes every
    # legend group's header AND item text (Sample Type, Peptide Counts,
    # Average Abundance) uniformly, in one shared Plotly legend.
    font_size_xaxis_label=14, font_size_yaxis_label=14,
    font_size_legend=14, font_size_plot_title=16,
    font_size_xaxis_tick=14, font_size_yaxis_tick=14,
    font_size_var_label=14,     # per-row peptide-count heatmap label (e.g. "High NE")
    legend_position='right',    # 'right' | 'below' — see LEGEND_POSITIONS
):
    """
    Plotly interactive version of the landscape / portrait heatmap.

    Layout (shared x-axis, top → bottom):
      Row 1       — Line-plot panel (averaged MS lines and/or individual peptide lines)
      Rows 2…n+1  — Heatmap-bar panel per variable
                      bar colour  = peptide count  (cmap)
                      bar height  = averaged abundance

    Hover text
      Line traces       : Position | Amino Acid | Abundance | Variable
                          [+ Peptide | Function when bioactive mode]
      Heatmap bar traces: Position | Amino Acid | Peptide Count | Abundance

    Returns a plotly.graph_objects.Figure, or None on failure.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        from matplotlib.colors import to_rgba
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("plotly is required for interactive heatmap: pip install plotly")
        return None

    var_keys = list(available_data_variables_dict.keys())
    num_vars = len(var_keys)
    if num_vars == 0:
        return None

    # ── y-axis tick values ────────────────────────────────────────────────────
    tick1 = float(y_ticks[0])  if y_ticks                else 0.0
    tick2 = float(y_ticks[1])  if len(y_ticks) > 1       else (tick1 + (float(y_ticks[-1]) if y_ticks else 1.0)) / 2
    tick3 = float(y_ticks[-1]) if y_ticks                else 1.0

    max_seq_len = max(
        len(available_data_variables_dict[v].get('AA_list', []))
        for v in var_keys
    )

    # ── decide whether to show the plain protein-sequence row ────────────────
    # Only shown when all variables share one protein.
    _prot_ids = list(dict.fromkeys(
        vd.get('protein_id') or vd.get('protein_name', '')
        for vd in available_data_variables_dict.values()
    ))
    _prot_ids = [p for p in _prot_ids if p]
    single_protein = len(_prot_ids) <= 1

    # protein AA list for the sequence strip (use first variable's AA_list)
    protein_aa_list = []
    if single_protein and var_keys:
        protein_aa_list = list(available_data_variables_dict[var_keys[0]].get('AA_list', []))

    # ── wrap width: the ONLY thing that separates Landscape from Portrait ─────
    # Mirrors the matplotlib renderers, where landscape uses
    # chunk_size = max_sequence_length (one band spanning the whole protein) and
    # portrait uses chunk_size = 78 (the sequence wrapped into stacked bands).
    # Same band-drawing code below serves both — Landscape is just the n_chunks==1
    # special case, so its output is unchanged from the pre-portrait version.
    PORTRAIT_CHUNK = 78
    is_portrait = str(orientation).strip().lower() == 'portrait'
    chunk_width = PORTRAIT_CHUNK if is_portrait else max_seq_len
    if chunk_width < 1:
        chunk_width = max_seq_len or 1
    n_chunks = max(1, -(-max_seq_len // chunk_width))            # ceil division
    # 1-indexed inclusive [start, end] position window per band
    chunk_windows = [
        (c * chunk_width + 1, min((c + 1) * chunk_width, max_seq_len))
        for c in range(n_chunks)
    ]

    # ── helper: matplotlib colour → Plotly rgba string ───────────────────────
    # Accepts colour names ('white'), hex, or RGBA tuples/arrays. The string
    # short-circuit must be isinstance-guarded: individual-peptide lines pass a
    # numpy RGBA array, and `array == 'white'` raises "truth value ambiguous"
    # (to_rgba handles both forms, so the guard just avoids that comparison).
    def _rgba(mpl_color):
        if isinstance(mpl_color, str) and mpl_color == 'white':
            return 'rgba(255,255,255,1.0)'
        try:
            r, g, b, a = to_rgba(mpl_color)
            return f'rgba({int(r*255)},{int(g*255)},{int(b*255)},{a:.3f})'
        except Exception:
            return 'rgba(180,180,180,1.0)'

    # ── colour map objects ────────────────────────────────────────────────────
    avg_cmap_obj = avg_cmap if not isinstance(avg_cmap, str) else plt.get_cmap(avg_cmap)
    lp_cmap_obj  = plt.get_cmap(lp_selected_color) if isinstance(lp_selected_color, str) else lp_selected_color

    # Legend group headers. Mirrors the matplotlib legend, which heads the
    # averaged / individual sample-type lines with legend_title[0] ("Sample Type:")
    # and the peptide-count swatches with legend_title[1]. The abundance scale
    # (legend_title[2]) is shown directly on the row-1 y-axis, not duplicated here.
    sample_type_title = legend_title[0] if legend_title else 'Sample Type:'

    # Legend placement ('right' = one combined legend beside the plot; 'below' =
    # one boxed unit per group under the x-axis, one per row). _legend_kwargs()
    # returns the per-trace kwargs; _legend_units is filled in draw order.
    legend_position = normalize_legend_position(legend_position)
    legend_below = legend_position == 'below'
    _legend_kwargs, _legend_units = _make_legend_placer(legend_position, font_size_legend)

    # ══════════════════════════════════════════════════════════════════════════
    # PRECOMPUTE — build every trace once at full sequence length, then slice it
    # per band in the render loop below. Keeping this position-indexed (rather
    # than re-deriving per band) is what lets one code path serve both layouts.
    # ══════════════════════════════════════════════════════════════════════════

    # line_specs: one entry per line trace (averaged or individual), each holding
    # its full-length x/y/hover plus styling. `is_indiv` only affects line width.
    line_specs = []   # dict(x, y, hover, name, color, legendgroup, is_indiv)

    # ── differential-comparison contrast line (replaces abundance lines) ─────
    # When a comparison track is supplied it is the ONLY row-1 trace: a single
    # signed line (positive => higher in Series A) over the density heatmap. The
    # abundance / individual-peptide lines are suppressed so the panel is
    # unambiguous. Undefined positions are NaN → Plotly renders them as gaps.
    contrast_is_active = comparison_track is not None
    if contrast_is_active:
        from data_analysis.services.stats import effect_size_label
        _track   = list(comparison_track.get('values', []))
        _metric  = comparison_track.get('metric', 'smd')
        _lbl_a   = comparison_track.get('label_a', 'A')
        _lbl_b   = comparison_track.get('label_b', 'B')
        _aa_list = list(comparison_track.get('aa_list', []))
        _x_vals  = list(range(1, len(_track) + 1))
        _y_vals  = [float(v) if v is not None and np.isfinite(v) else None for v in _track]
        # The metric is already named on the y-axis, so the legend only needs to
        # say WHAT is being compared ("Threshold vs Low") under a "Comparison"
        # header. Hover keeps a short metric tag for per-residue readout.
        _metric_short = 'SMD' if _metric == 'smd' else 'log₂FC'
        _trace_name   = f'{_lbl_a} vs {_lbl_b}'
        _hover = []
        for _x, _y in zip(_x_vals, _y_vals):
            _aa = _aa_list[_x - 1] if 0 < _x <= len(_aa_list) else '?'
            if _y is None:
                _hover.append(f"Position: {_x}<br>Amino Acid: {_aa}<br>{_metric_short}: n/a")
            elif _metric == 'smd':
                _hover.append(
                    f"Position: {_x}<br>Amino Acid: {_aa}<br>"
                    f"{_metric_short}: {_y:+.2f} ({effect_size_label(_y)})<br>"
                    f"(+ = higher in {_lbl_a})")
            else:
                _hover.append(
                    f"Position: {_x}<br>Amino Acid: {_aa}<br>"
                    f"{_metric_short}: {_y:+.2f}<br>(+ = higher in {_lbl_a})")
        line_specs.append(dict(
            x=_x_vals, y=_y_vals, hover=_hover, name=_trace_name,
            color='rgba(33,33,33,0.9)', legendgroup='contrast', is_indiv=False,
            grouptitle='Comparison',
        ))

    # ── averaged MS lines (one per variable) ─────────────────────────────────
    if not contrast_is_active and ms_average_choice in ('yes', 'only'):
        for var_idx, var_key in enumerate(var_keys):
            vd      = available_data_variables_dict[var_key]
            ms_data = vd.get('ms_data_list', [])
            aa_list = vd.get('AA_list', [])
            label   = vd.get('label', var_key)
            color   = _rgba(avg_cmap_obj(var_idx % avg_cmap_obj.N))
            if isinstance(ms_data, (pd.DataFrame, pd.Series)):
                x_vals = list(ms_data.index)
                y_vals = list(ms_data.values)
            else:
                x_vals = list(range(1, len(ms_data) + 1))
                y_vals = list(ms_data)

            aa_at = [aa_list[x - 1] if 0 < x <= len(aa_list) else '?' for x in x_vals]
            try:
                hover = [
                    f"Position: {x}<br>Amino Acid: {aa}<br>Variable: {label}<br>Abundance: {float(y):.3e}"
                    for x, aa, y in zip(x_vals, aa_at, y_vals)
                ]
            except Exception:
                hover = [f"Position: {x}" for x in x_vals]
            line_specs.append(dict(
                x=x_vals, y=y_vals, hover=hover, name=label,
                color=color, legendgroup='avg_lines', is_indiv=False,
            ))

    # ── individual peptide / bioactive-function lines ─────────────────────────
    if not contrast_is_active \
            and (bio_or_pep != 'no' or filter_type in ('bioactive-only', 'functional-only')) \
            and ms_average_choice != 'only':
        for var_idx, var_key in enumerate(var_keys):
            vd      = available_data_variables_dict[var_key]
            label   = vd.get('label', var_key)
            bp_abs  = vd.get('bioactive_peptide_abs_df')
            bp_func = vd.get('bioactive_peptide_func_df')
            aa_list = vd.get('AA_list', [])

            if bp_abs is None or bp_abs.empty:
                continue

            # apply the same filter logic as the Matplotlib functions
            filtered_abs, filtered_func, _ = filter_data_by_selection(
                bp_abs, bp_func, selected_peptides, selected_functions,
                bio_or_pep, filter_type, ms_average_choice, label,
            )
            if filtered_abs is None or filtered_abs.empty:
                continue

            # resolve effective function list (mirrors landscape logic)
            eff_bio    = bio_or_pep
            eff_funcs  = list(selected_functions)
            if filter_type == 'bioactive-only' and bio_or_pep == 'no':
                all_funcs = set()
                if filtered_func is not None and not filtered_func.empty:
                    for col in filtered_func.columns:
                        for val in filtered_func[col].dropna():
                            for part in str(val).split(';'):
                                part = part.strip()
                                if part and part != 'nan':
                                    all_funcs.add(part)
                if all_funcs:
                    eff_bio   = '2'
                    eff_funcs = sorted(all_funcs)

            # colour map for individual lines
            items = (list(selected_peptides) if eff_bio == '1' else eff_funcs) if eff_bio != 'no' else []
            n_items    = max(len(items), 1)
            lp_colors  = lp_cmap_obj(np.linspace(0, 1, n_items))
            item_cmap  = {item: lp_colors[i % n_items] for i, item in enumerate(items)}
            item_cmap['Multiple'] = (0.0, 0.0, 0.0, 1.0)

            for col in filtered_abs.columns:
                if col == 'average':
                    continue
                y_series = numeric_column(filtered_abs, col).replace(0, np.nan)
                if y_series.isna().all():
                    continue

                x_vals = list(y_series.index)
                y_vals = list(y_series.values)
                aa_at  = [aa_list[x - 1] if 0 < x <= len(aa_list) else '?' for x in x_vals]

                # determine label and function string
                pep_label = col
                func_str  = ''
                if eff_bio == '2' and filtered_func is not None and col in filtered_func.columns:
                    fvals = filtered_func[col].dropna()
                    if not fvals.empty:
                        found = set()
                        for val in fvals:
                            for sub in str(val).split(';'):
                                sub = sub.strip()
                                for func in eff_funcs:
                                    if func.lower().replace(' ', '') == sub.lower().replace(' ', ''):
                                        found.add(func)
                        if len(found) > 1:
                            pep_label = 'Multiple'
                            func_str  = '; '.join(sorted(found))
                        elif found:
                            pep_label = next(iter(found))
                            func_str  = pep_label

                line_color = _rgba(item_cmap.get(pep_label, (0.5, 0.5, 0.5, 1.0)))
                # Hover kept aligned 1:1 with x/y (no NaN filtering) so a band
                # can slice x/y/hover together by index without drift.
                hover = [
                    f"Position: {x}<br>Amino Acid: {aa}<br>Variable: {label}"
                    f"<br>Abundance: {y:.3e}<br>Peptide: {col}"
                    + (f"<br>Function: {func_str}" if func_str else "")
                    for x, aa, y in zip(x_vals, aa_at, y_vals)
                ]
                line_specs.append(dict(
                    x=x_vals, y=y_vals, hover=hover, name=pep_label,
                    color=line_color, legendgroup='pep_lines', is_indiv=True,
                ))

    # ── plain protein-sequence strip (single protein only) ───────────────────
    seq_positions = []
    seq_text_full = []
    seq_hover_full = []
    if single_protein and protein_aa_list:
        seq_positions  = list(range(1, len(protein_aa_list) + 1))
        seq_text_full  = list(protein_aa_list)
        seq_hover_full = [f"Position: {p}<br>Amino Acid: {aa}"
                          for p, aa in zip(seq_positions, protein_aa_list)]

    # ── colored heatmap tile rows (one spec per variable) ─────────────────────
    hm_specs = []   # dict(var_key, label, positions, bar_colors, hover)
    for i, var_key in enumerate(var_keys):
        vd          = available_data_variables_dict[var_key]
        hm_label    = vd.get('label', var_key)

        # --- safe extraction: avoid `series or []` which raises ambiguous-truth error ---
        _raw_counts = vd.get('peptide_counts')
        if _raw_counts is None:
            counts_list = []
        else:
            try:
                counts_list = list(_raw_counts)
            except Exception:
                counts_list = []

        _raw_ms = vd.get('ms_data_list')
        ms_data = _raw_ms if _raw_ms is not None else []

        aa_list = vd.get('AA_list', [])

        if len(counts_list) == 0:
            hm_df = vd.get('heatmap_df')
            if hm_df is not None and not hm_df.empty:
                counts_list = hm_df['count'].tolist()
                ms_data     = hm_df['average'].fillna(0).tolist()
                aa_list     = hm_df['AA'].tolist()

        if len(counts_list) == 0:
            continue

        positions      = list(range(1, len(counts_list) + 1))
        bar_colors_mpl = get_grouped_colors(counts_list, max_count, num_colors, plot_zero, cmap)
        bar_colors     = [_rgba(c) for c in bar_colors_mpl]

        # --- normalise ms_data to a plain Python list ---
        if isinstance(ms_data, (pd.DataFrame, pd.Series)):
            ms_data = list(ms_data.values if hasattr(ms_data, 'values') else ms_data)
        elif not isinstance(ms_data, list):
            try:
                ms_data = list(ms_data)
            except Exception:
                ms_data = [0.0] * len(counts_list)

        # safe float helper for abundance hover text
        def _safe_float(v):
            try:
                f = float(v)
                return f if f == f else 0.0   # NaN → 0
            except (TypeError, ValueError):
                return 0.0

        # abundance values for hover (ms_data may be shorter than counts if fallback used)
        abun_vals = [_safe_float(ms_data[k]) if k < len(ms_data) else 0.0
                     for k in range(len(counts_list))]

        # ── Build per-position peptide overlap map ────────────────────────────
        # Maps 1-indexed position → list of (peptide_seq, start_pos_1indexed).
        # Used to render the aligned sequence block in hover text.
        _pep_hover_map = {}
        _bp_abs_h  = vd.get('bioactive_peptide_abs_df')
        _bp_func_h = vd.get('bioactive_peptide_func_df')
        if _bp_abs_h is not None and not _bp_abs_h.empty:
            try:
                _fabs, _, _ = filter_data_by_selection(
                    _bp_abs_h, _bp_func_h,
                    selected_peptides, selected_functions,
                    bio_or_pep, filter_type, ms_average_choice, hm_label,
                )
                if _fabs is None or _fabs.empty:
                    _fabs = _bp_abs_h
            except Exception:
                _fabs = _bp_abs_h
            for _col in list(_fabs.columns):
                if _col in ('average', 'count', 'AA'):
                    continue
                _ys = numeric_column(_fabs, _col).replace(0, np.nan)
                _valid = _ys.notna()
                if not _valid.any():
                    continue
                # DataFrame index is 0-based; positions used in hover text are
                # 1-based (range(1, …)).  Convert here so every key and _sp
                # value in the map is 1-indexed, matching the pos argument
                # received by _pep_block.
                _sp = int(_valid[_valid].index.min()) + 1   # → 1-indexed start
                for _p in _valid[_valid].index:
                    _pep_hover_map.setdefault(int(_p) + 1, []).append((_col, _sp))

        _MAX_PEP = 25   # cap to keep hover text manageable

        def _pep_block(pos, _pep_hover_map=_pep_hover_map, _MAX_PEP=_MAX_PEP):
            raw = _pep_hover_map.get(int(pos), [])
            # Keep only peptides where the hovered position genuinely falls within
            # the peptide.  Column names have the form "start-end" or "start-end UID";
            # parse the true peptide length (end - start + 1) from that prefix rather
            # than using len(col_label), which is the string length of the label and
            # is almost always shorter than the actual residue span — causing the last
            # several positions of every long peptide to be silently dropped.
            def _pep_len(col):
                m = re.match(r'^(\d+)-(\d+)', col)
                return int(m.group(2)) - int(m.group(1)) + 1 if m else len(col)

            entries = [(seq, st) for seq, st in raw
                       if 0 <= int(pos) - st < _pep_len(seq)]
            if not entries:
                return ''
            orig_len   = len(entries)
            truncated  = orig_len > _MAX_PEP
            shown      = entries[:_MAX_PEP]
            offsets    = [int(pos) - s for _, s in shown]
            max_off    = max(offsets)
            n_dig      = len(str(orig_len))   # for consistent ID width
            inner = []
            for uid, ((seq, _st), off) in enumerate(zip(shown, offsets), start=1):
                pad = max_off - off
                nbsp = '&nbsp;' * pad
                # Column label format: "start-end" or "start-end UID_SEQUENCE".
                # The amino acid at peptide offset `off` sits inside the UID part,
                # which begins after the "start-end " prefix.  Compute that prefix
                # length so we bold the correct character rather than a digit/dash.
                _pfx_m = re.match(r'^\d+-\d+\s+', seq)
                _pfx   = len(_pfx_m.group(0)) if _pfx_m else len(seq)
                _char_idx = _pfx + off
                if 0 <= _char_idx < len(seq):
                    # Plotly hover supports only <b><i><br> etc — <u> renders
                    # as literal text.  Use square brackets + bold instead.
                    seq_html = (nbsp + seq[:_char_idx]
                                + f'[<b>{seq[_char_idx]}</b>]'
                                + seq[_char_idx + 1:])
                else:
                    seq_html = nbsp + seq   # no sequence in label — show as-is
                id_str = str(uid).rjust(n_dig)
                inner.append(f'{id_str}.&nbsp;{seq_html}')
            if truncated:
                inner.append(f'&nbsp;&nbsp;&nbsp;… and {orig_len - _MAX_PEP} more')
            block = '<br>'.join(inner)
            return (
                '<br><b>Overlapping Peptides:</b><br>'
                f'<span style="font-family:monospace">{block}</span>'
            )

        # AA letters are already printed in the plain-sequence strip above, so we
        # do NOT repeat them on the colored tile rows. Hover carries the detail.
        aa_full = [aa_list[k] if k < len(aa_list) else '?' for k in range(len(counts_list))]
        hover = [
            f"<b>{hm_label}</b><br>"
            f"<b>Position: {pos} | {aa}</b><br>"
            f"Peptide Count: {cnt}&nbsp;&nbsp;Averaged Abundance: {abun:.3e}"
            + _pep_block(pos)
            for pos, aa, cnt, abun in zip(positions, aa_full, counts_list, abun_vals)
        ]

        hm_specs.append(dict(
            var_key=var_key, label=hm_label,
            positions=positions, bar_colors=bar_colors, hover=hover,
            aa=aa_full,
        ))

    # ══════════════════════════════════════════════════════════════════════════
    # ROW PLAN — one band per chunk; within a band the row order matches the
    # matplotlib layout (line-plot → AA-sequence strip → one tile row per var).
    # Landscape has a single band so this collapses to the original layout.
    # ══════════════════════════════════════════════════════════════════════════
    LINE_PX = 260       # line-plot row height
    SEQ_PX  = 25        # plain AA-sequence strip
    HM_PX   = 26        # each colored heatmap tile row
    # Small inter-band gap: just enough to clear each band's position ladder so
    # the wrapped bands still read as one continuous figure (was 46).
    GAP_PX  = 26

    # The shared grey AA-sequence strip is shown only for a single protein AND
    # only when the letters are NOT being printed on the tiles themselves. With
    # aa_on_tiles the sequence lives on each sample's row instead (and multi-protein
    # comparisons — which have no shared strip — get a per-protein sequence too).
    show_seq = bool(single_protein and seq_positions and not aa_on_tiles)

    rows_meta = []      # (kind, chunk_idx, var_key_or_None)
    for c in range(n_chunks):
        if c > 0:
            rows_meta.append(('gap', c, None))
        rows_meta.append(('line', c, None))
        if show_seq:
            rows_meta.append(('seq', c, None))
        for spec in hm_specs:
            rows_meta.append(('hm', c, spec['var_key']))

    _ROW_PX = {'line': LINE_PX, 'seq': SEQ_PX, 'hm': HM_PX, 'gap': GAP_PX}
    num_rows   = len(rows_meta)
    px_heights = [_ROW_PX[r[0]] for r in rows_meta]
    total_px   = sum(px_heights)
    row_heights_frac = [h / total_px for h in px_heights]

    # paper-space domain centres for left-side annotations (rows go top → bottom)
    _cumfrac   = [sum(row_heights_frac[:k]) for k in range(num_rows + 1)]
    row_y_top    = [1.0 - _cumfrac[k]     for k in range(num_rows)]
    row_y_bottom = [1.0 - _cumfrac[k + 1] for k in range(num_rows)]
    row_y_center = [(row_y_top[k] + row_y_bottom[k]) / 2 for k in range(num_rows)]

    def _row_num(kind, chunk, var_key=None):
        for k, (rk, rc, rv) in enumerate(rows_meta):
            if rk == kind and rc == chunk and (var_key is None or rv == var_key):
                return k + 1        # 1-indexed for Plotly
        return None

    fig = make_subplots(
        rows=num_rows, cols=1,
        shared_xaxes=False,     # per-band x ranges (each chunk shows its own window)
        vertical_spacing=0,     # zero gap → rows touch; bands separated by 'gap' rows
        row_heights=row_heights_frac,
    )

    # ── abundance log transform ──────────────────────────────────────────────
    # When the user asks for a log transform, transform the plotted abundance
    # DATA (line heights) — not just relabel a log-scaled axis — so the bars /
    # lines scale with the transform and log10 vs log2 are visibly different.
    # The axis is then plain linear over the transformed values and the ticks
    # read as log10 / log2 values, never raw abundance.
    _do_log = bool(log_transform) and not contrast_is_active and tick1 > 0
    _log_sub_str = _log_sub(log_base)

    if _do_log:
        for spec in line_specs:
            new_y, new_h = [], []
            for yi, hi in zip(spec['y'], spec['hover']):
                ly = _log_value(yi, log_base)
                new_y.append(ly)
                if isinstance(hi, str):
                    tag = (f"<br>Log<sub>{_log_sub_str}</sub> Abundance: {ly:.2f}"
                           if ly is not None else
                           f"<br>Log<sub>{_log_sub_str}</sub> Abundance: n/a")
                    new_h.append(hi + tag)
                else:
                    new_h.append(hi)
            spec['y'] = new_y
            spec['hover'] = new_h
        _lt1, _lt2, _lt3 = (_log_value(tick1, log_base),
                            _log_value(tick2, log_base),
                            _log_value(tick3, log_base))

    # y-axis tick label formatter (3 sig-fig scientific for raw abundance; the
    # transformed-value case is formatted inline below).
    def _fmt_tick(v):
        if v == 0:
            return '0'
        return f"{v:.3g}"

    if _do_log:
        lp_range     = [_lt1, _lt3]
        lp_type      = 'linear'
        lp_tick_vals = [_lt1, _lt2, _lt3]
        lp_tick_text = [f"{v:.2f}" if v is not None else '' for v in lp_tick_vals]
    else:
        lp_range     = [tick1, tick3]
        lp_type      = 'linear'
        lp_tick_vals = [tick1, tick2, tick3]
        lp_tick_text = [_fmt_tick(tick1), _fmt_tick(tick2), _fmt_tick(tick3)]

    x_label_str      = xaxis_label or 'Protein Sequence Position'
    legend_names_shown = set()   # show each legend entry once, figure-wide

    # ── compact left-side labelling for tall wrapped portraits ────────────────
    # When a long protein wraps into many stacked bands, repeating the
    # per-variable ("High NE" / "Low NE") row label on every band packs them
    # so tightly they overrun into an illegible column. At 1-4 bands the
    # repetition is harmless (arguably clearer, if redundant) and stays on
    # every band; from 5 bands on it's dropped down to ONE label (on the top
    # band), letting position alone carry the identity of every band below.
    # Landscape (always n_chunks == 1) always keeps its per-band label.
    compact_labels = is_portrait and n_chunks > 4

    # Per-band x-axis tick-label font: a fixed size reads fine for a handful
    # of bands, but a long protein wraps into dozens of stacked bands, each
    # repeating its own position ladder — shrink the font as bands accumulate
    # so the figure doesn't get visually noisy. Floors at 8px. Landscape is
    # always a single band, so it keeps the base size.
    _x_tick_font_size = max(8, font_size_xaxis_tick - (n_chunks - 1) // 3) if is_portrait else font_size_xaxis_tick

    # Peptide-count row ("High NE" / "Low NE") label font: same reasoning —
    # at 1-4 bands compact_labels is off so the label repeats on EVERY band,
    # and each repeat is drawn at a fixed pixel distance from the plot edge
    # regardless of band count, so more bands means more of these labels
    # packed down the same-width left margin. Shrinks off the user's base size
    # as bands accumulate, same as the x-axis tick font above.
    _var_label_font_size = max(8, font_size_var_label - (n_chunks - 1) // 3) if is_portrait else font_size_var_label

    # ══════════════════════════════════════════════════════════════════════════
    # RENDER — draw every band by slicing the precomputed traces to its window.
    # ══════════════════════════════════════════════════════════════════════════
    for c, (win_start, win_end) in enumerate(chunk_windows):
        line_row = _row_num('line', c)
        band_rows = [k + 1 for k, (rk, rc, rv) in enumerate(rows_meta)
                     if rc == c and rk != 'gap']
        band_bottom = band_rows[-1]     # last hm row of the band → carries x ticks

        # ── line-plot traces (sliced to the window) ──────────────────────────
        for spec in line_specs:
            xs, ys, hs = [], [], []
            for xi, yi, hi in zip(spec['x'], spec['y'], spec['hover']):
                try:
                    xf = float(xi)
                except (TypeError, ValueError):
                    continue
                if win_start - 0.5 <= xf <= win_end + 0.5:
                    xs.append(xi); ys.append(yi); hs.append(hi)
            if not xs:
                continue
            show_leg = spec['name'] not in legend_names_shown
            if show_leg:
                legend_names_shown.add(spec['name'])
            fig.add_trace(
                go.Scatter(
                    x=xs, y=ys, mode='lines', name=spec['name'],
                    line=dict(color=spec['color'], width=1.2 if spec['is_indiv'] else 1.5),
                    hoverinfo='text', hovertext=hs,
                    **_legend_kwargs(spec.get('grouptitle', sample_type_title),
                                     spec['legendgroup'],
                                     entry=spec['name'] if show_leg else None),
                    showlegend=show_leg,
                ),
                row=line_row, col=1,
            )

        # line-plot y-axis (abundance) — one per band. The axis TITLE is drawn
        # once as a single figure-level label (below) rather than repeated on
        # every band. Portrait mirrors the Compact renderer: numeric tick
        # labels are dropped from EVERY band (they'd repeat down dozens of
        # stacked bands for a long protein) and shown once instead as
        # Min/Mid/Max entries in the legend — see below. Every band still
        # carries gridlines at the same three values plus the zero reference
        # in comparison mode. Landscape (single band) keeps its on-axis ticks.
        _show_yticklabels = not is_portrait
        fig.update_yaxes(
            title=None,
            type=lp_type, range=lp_range,
            # Zoom is sequence-direction only: the abundance scale is shared by
            # every band and doubles as the legend's Min/Mid/Max key, so letting
            # the user rescale it (modebar zoom-in / box-zoom / scroll) would
            # silently break that correspondence and clip the trace. See the
            # matching fixedrange on the tile rows below.
            fixedrange=True,
            showgrid=True, gridwidth=1, gridcolor='rgba(200,200,200,0.5)',
            tickvals=lp_tick_vals,
            ticktext=lp_tick_text,
            showticklabels=_show_yticklabels,
            ticks='outside' if _show_yticklabels else '',
            ticklen=6, tickwidth=2,
            tickfont=dict(size=font_size_yaxis_tick),
            # Signed differential track: draw a solid zero reference line so
            # positive (higher in A) / negative (higher in B) excursions read
            # against a fixed baseline. Positive-only abundance path unchanged.
            zeroline=contrast_is_active,
            zerolinecolor='rgba(80,80,80,0.9)' if contrast_is_active else None,
            zerolinewidth=1.5 if contrast_is_active else None,
            row=line_row, col=1,
        )

        # ── plain protein-sequence strip (sliced) ────────────────────────────
        if show_seq:
            seq_row = _row_num('seq', c)
            i0, i1 = win_start - 1, win_end          # contiguous slice
            pos_s  = seq_positions[i0:i1]
            aa_s   = seq_text_full[i0:i1]
            hov_s  = seq_hover_full[i0:i1]
            show_seq_text = len(pos_s) <= 300
            text_s = aa_s if show_seq_text else [''] * len(pos_s)
            fig.add_trace(
                go.Bar(
                    x=pos_s, y=[1.0] * len(pos_s),
                    text=text_s, textposition='inside',
                    # AA-lettering size is per-residue-cell constrained, not a free
                    # legibility knob: Landscape puts the whole protein (up to 200+
                    # residues) in one band, so px-per-residue is already tight at
                    # size 13 — going bigger overlaps adjacent letters. Portrait
                    # always wraps to a fixed PORTRAIT_CHUNK=78-residue band
                    # regardless of protein length, so for anything longer than 78
                    # residues it gets proportionally more px per residue and can
                    # afford a larger letter. Ed-14: raising this declared size
                    # alone does NOT grow the rendered glyphs — Plotly's
                    # uniformtext shrinks text to fit the per-residue cell width
                    # regardless of the declared size, and cell width (driven by
                    # PORTRAIT_CHUNK) is already the binding constraint at 16px.
                    # Confirmed empirically: bumping to 26px produced pixel-
                    # identical output. Left at 16px; the real lever is chunk width.
                    textfont=dict(size=(16 if is_portrait else 13), color='#333333'),
                    insidetextanchor='middle', textangle=0,
                    marker_color='rgba(230,230,230,1.0)', marker_line_width=0,
                    width=1.0, hoverinfo='text', hovertext=hov_s,
                    showlegend=False, name='Protein Sequence',
                ),
                row=seq_row, col=1,
            )
            fig.update_yaxes(range=[0, 1], showgrid=False, showticklabels=False,
                             ticks='', fixedrange=True, row=seq_row, col=1)

        # ── colored heatmap tile rows (one per variable, sliced) ─────────────
        for spec in hm_specs:
            hm_row_num = _row_num('hm', c, spec['var_key'])
            i0, i1 = win_start - 1, win_end
            pos_s   = spec['positions'][i0:i1]
            col_s   = spec['bar_colors'][i0:i1]
            hov_s   = spec['hover'][i0:i1]
            if pos_s:
                # With aa_on_tiles, print the AA letter on each tile (so each
                # sample carries its own sequence — ideal for comparing proteins
                # or variants); otherwise the tiles are blank and the shared grey
                # strip above carries the sequence.
                tile_text = spec['aa'][i0:i1] if aa_on_tiles else [''] * len(pos_s)
                fig.add_trace(
                    go.Bar(
                        x=pos_s, y=[1.0] * len(pos_s),
                        text=tile_text, textposition='inside',
                        # See the sequence-strip trace above: same Landscape-tight /
                        # Portrait-has-room reasoning applies to on-tile AA letters.
                        textfont=dict(size=(15 if is_portrait else 12), color='black'),
                        insidetextanchor='middle', textangle=0,
                        marker_color=col_s, marker_line_width=0, width=1.0,
                        hoverinfo='text', hovertext=hov_s,
                        showlegend=False, name=spec['label'],
                    ),
                    row=hm_row_num, col=1,
                )
            # [0, 1] is a pure layout device (every tile is a full-height bar) —
            # there is nothing to zoom into vertically, and a y-zoom here just
            # slices the tiles in half.
            fig.update_yaxes(type='linear', range=[0, 1], showgrid=False,
                             showticklabels=False, ticks='', fixedrange=True,
                             row=hm_row_num, col=1)
            # variable-name label on the left (mirrors ax.text(..., ha='right')).
            # Skip when this variable contributes no residues to the band (ragged
            # multi-protein: a shorter sequence doesn't reach later chunks) so the
            # label never floats over an empty row. In compact mode the variable
            # is named once on the top band (position identifies the rest).
            if pos_s and ((not compact_labels) or (c == 0)):
                # Pinned 6px left of the plot edge in PIXEL space (xshift) so the
                # label sits at the same place regardless of browser width.
                fig.add_annotation(
                    text=spec['label'], xref='paper', yref='paper',
                    x=0, y=row_y_center[hm_row_num - 1], xshift=-6,
                    xanchor='right', yanchor='middle',
                    showarrow=False, font=dict(size=_var_label_font_size),
                )

        # ── per-band x-axis: absolute-position ladder on the band's bottom row ─
        # Pad the left so the first bar (centre = win_start, left edge −0.5) is
        # fully visible; clamp pan/zoom to the band window (Plotly ≥ 5.17).
        # Universal scale: every band spans a FIXED chunk_width residues (not
        # win_end), so bars and AA letters render at the same size across all
        # bands. A short final band simply leaves empty space on the right rather
        # than stretching its few residues across the full width.
        x_min = win_start - 2.0
        x_max = win_start + chunk_width - 1 + 0.5
        x_common = dict(
            autorange=False, range=[x_min, x_max],
            minallowed=x_min, maxallowed=x_max,
            showgrid=True, gridwidth=1, gridcolor='rgba(200,200,200,0.4)',
            zeroline=False,
        )
        # Every row in the band is a window onto the SAME residue axis (line
        # plot, AA strip, tile rows), so they must zoom and pan as one unit —
        # otherwise zooming the abundance trace leaves the tiles beneath it at
        # the old scale and the figure stops lining up. `matches` makes the
        # band's bottom row (the one carrying the tick ladder) the master.
        # Bands are NOT matched to each other: each shows a different residue
        # window, so a shared range would be meaningless.
        _x_master = 'x' if band_bottom == 1 else f'x{band_bottom}'
        for r in band_rows:
            if r == band_bottom:
                # x-axis TITLE only on the final band; every band shows its ladder.
                title = dict(text=x_label_str, font=dict(size=font_size_xaxis_label)) if c == n_chunks - 1 else None
                fig.update_xaxes(
                    title=title, dtick=20, tick0=0, tickfont=dict(size=_x_tick_font_size),
                    showticklabels=True, **x_common, row=r, col=1,
                )
            else:
                fig.update_xaxes(showticklabels=False, matches=_x_master,
                                 **x_common, row=r, col=1)

        # ── silence the spacer 'gap' row above this band (portrait only) ──────
        if c > 0:
            gap_row = _row_num('gap', c)
            fig.update_xaxes(visible=False, row=gap_row, col=1)
            fig.update_yaxes(visible=False, range=[0, 1], showticklabels=False,
                             fixedrange=True, row=gap_row, col=1)

    # ── dynamic left margin + single shared y-axis title ─────────────────────
    # Paper-fraction offsets scale with plot width (tiny when the window is narrow
    # → the title overlaps the ticks / sample names; huge when full-screen → the
    # title clips off the left edge). So size the left margin to the longest
    # left-side label and pin the title a FIXED pixel distance from the figure's
    # left edge via `xshift`, making it width-invariant.
    _max_label_len = max((len(str(s['label'])) for s in hm_specs), default=0)
    # px-per-char scales with the actual label font (7.2 is the measured
    # width at the base size 12; _var_label_font_size shrinks it for tall
    # wrapped portraits, see above) so the reserved margin tracks the font.
    _var_label_px  = _max_label_len * (7.2 * _var_label_font_size / 12) + 6
    # Abundance ticks e.g. "1.8e+08" — Landscape only; Portrait drops them
    # from the axis entirely (moved to the legend, see below) so it doesn't
    # need this reserved.
    _tick_label_px = 0 if is_portrait else 54
    _label_zone_px = max(_var_label_px, _tick_label_px)
    _ytitle_px     = 28                            # rotated title band width (font 15)
    left_margin    = int(min(300, max(90, _ytitle_px + 14 + _label_zone_px)))

    # A user-supplied Y-Axis Label always wins. Otherwise fall back to the
    # effect-size metric name in comparison mode, or the abundance default.
    if yaxis_label:
        _ytitle = yaxis_label
    elif contrast_is_active:
        _ytitle = ('Standardized Mean Difference (Cohen’s d)'
                   if comparison_track.get('metric') == 'smd'
                   else 'log₂ Fold Change')
    elif _do_log:
        _ytitle = f'Log<sub>{_log_sub_str}</sub> Averaged Peptide Abundance'
    else:
        _ytitle = 'Averaged Peptide Abundance'
    fig.add_annotation(
        text=_ytitle,
        xref='paper', yref='paper', x=0, y=0.5,
        xshift=-(left_margin - 12), textangle=-90,
        xanchor='center', yanchor='middle',
        showarrow=False, font=dict(size=font_size_yaxis_label),
    )

    # ── legend: peptide-count colour scale only (abundance is on the y-axis) ──
    # Header rendered as a native Plotly legend-group title (like the Sample Type
    # group above), not a bold placeholder legend item — so it aligns flush-left
    # (no marker-column indent) and shares one grouptitlefont with every header.
    legend_title_pep = legend_title[1] if len(legend_title) > 1 else 'Peptide Counts:'

    hm_handles, hm_labels = create_heatmap_legend_handles(cmap, num_colors, max_count, plot_zero)
    for lbl, handle in zip(hm_labels, hm_handles):
        fc = handle.get_facecolor()
        r_c, g_c, b_c, a_c = fc
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(
                size=14,
                color=f'rgba({int(r_c*255)},{int(g_c*255)},{int(b_c*255)},{a_c:.3f})',
                symbol='square', line=dict(width=1, color='black'),
            ),
            name=lbl, showlegend=True,
            **_legend_kwargs(legend_title_pep, 'pep_count', entry=lbl),
        ))

    # ── legend: abundance/contrast scale reference (Portrait only) ───────────
    # Landscape shows Min/Mid/Max directly as numeric y-axis ticks (a single
    # band, so no repetition concern). Portrait can wrap a long protein into
    # dozens of bands, so — mirroring the Compact renderer's legend-based
    # abundance scale — those per-band numeric ticks are dropped (see the
    # y-axis update above) and shown once here instead.
    if is_portrait:
        _abun_legend_title = (
            f"{_ytitle}:" if contrast_is_active
            else (legend_title[2] if len(legend_title) > 2 else 'Average Abundance:')
        )
        if _do_log and not contrast_is_active:
            _abun_legend_title = f"Log<sub>{_log_sub_str}</sub> {_abun_legend_title}"
        for tick_disp, tick_lbl in zip(lp_tick_text, ['Min', 'Mid', 'Max']):
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(size=10, symbol='line-ew', color='black', line=dict(width=2)),
                name=f"{tick_lbl}: {tick_disp}",
                showlegend=True,
                **_legend_kwargs(_abun_legend_title, 'abundance_ticks',
                                 entry=f"{tick_lbl}: {tick_disp}"),
            ))

    # ── overall layout ────────────────────────────────────────────────────────
    # No figure title: the protein is already named by the "<name> Sequence"
    # x-axis and the page heading, so a plot title only repeats it. The original
    # notebook heatmaps carried none either.

    # Optional user title (blank → Plotly draws none). Reserve top margin only
    # when a title is present so the untitled default stays tight.
    title_str = (plot_title or '').strip()

    # ── margins + legend geometry ────────────────────────────────────────────
    # The 'below' modes free the right margin entirely (that width goes to the
    # heatmap) and buy the reclaimed space back at the bottom for the legend
    # strip, which grows with the number of rows the units pack into.
    top_margin   = 60 if title_str else 30
    right_margin = 40 if legend_below else 180

    # figure height with NO legend strip = every row's pixel height + margins
    base_height  = total_px + top_margin + LEGEND_BASE_BOTTOM_PX + 10
    legend_meta  = None

    # Portrait wraps the sequence into many stacked bands, so a comparison's
    # zero-line ticks and/or a below-axis legend strip sit closer to the last
    # band than in Landscape (one band, more headroom). A little extra bottom
    # margin keeps them from crowding the final band. Added post-hoc to both
    # fig_height and bottom_margin so the plotted rows themselves are unchanged
    # — only the blank space below them grows.
    _extra_bottom_px = 20 if (is_portrait and (contrast_is_active or legend_below)) else 0

    if legend_below:
        legend_layout, fig_height, bottom_margin, legend_meta = _below_legend_layout(
            _below_legend_plan(_legend_units, font_size_legend),
            True, top_margin, base_height, font_size_legend)
    else:
        fig_height    = base_height
        bottom_margin = LEGEND_BASE_BOTTOM_PX

    fig_height    += _extra_bottom_px
    bottom_margin += _extra_bottom_px

    if not legend_below:
        legend_layout = dict(legend=dict(
            yanchor="top", y=0.99,
            xanchor="left", x=1.01,
            bordercolor="black", borderwidth=1,
            orientation="v",
            bgcolor="rgba(255,255,255,0.9)",
            # One shared font size for every legend group's header AND items —
            # Sample Type, Peptide Counts, and (Portrait only) Average Abundance
            # all live in this one Plotly legend.
            font=dict(size=font_size_legend),
            grouptitlefont=dict(size=font_size_legend, color='black'),
            itemsizing='constant',
            groupclick="toggleitem",
        ))

    if legend_meta is not None:
        # Lets the browser and the static-export view re-pack the units once the
        # real plotting width is known. See repack_below_legends().
        fig.update_layout(meta=legend_meta)

    fig.update_layout(
        title=dict(text=title_str, font=dict(size=font_size_plot_title, color='black'), x=0.5),
        # Explicit default size (not just family) so any text left unset above
        # (e.g. the row-1 abundance y-axis ticks) has a deterministic size
        # instead of silently inheriting Plotly's built-in default.
        font=dict(family=PLOTLY_FONT_FAMILY, color='black', size=12),
        height=fig_height,
        autosize=True,
        bargap=0, bargroupgap=0,
        plot_bgcolor='white',
        paper_bgcolor='white',
        showlegend=True,
        # left sized to longest sample name + y-title
        margin=dict(l=left_margin, r=right_margin, t=top_margin, b=bottom_margin),
        # Force bar text to always render, even when bars are narrower than the text.
        # Without this Plotly silently hides inside-text on narrow bars (e.g. 200+ AA sequence).
        uniformtext=dict(mode='show', minsize=10),
        **legend_layout,
    )
    return fig

"""_______________________________Compact Interactive Plot (Plotly)__________________________"""
def visualize_sequence_heatmap_compact(
    available_data_variables_dict,
    xaxis_label,
    yaxis_label,
    legend_title,
    cmap,
    y_ticks,
    log_transform,
    log_base=10,          # base for the abundance log transform (10 or 2)
    plot_title='',        # optional user figure title; blank → no title
    font_size_xaxis_label=14, font_size_yaxis_label=14,
    font_size_legend=14,
    font_size_plot_title=16, font_size_xaxis_tick=14,
    legend_position='right',   # 'right' | 'below' — see LEGEND_POSITIONS
):
    """
    Generate an interactive Plotly compact heatmap.

    Bar height  = averaged peptide abundance (ms_data_list).
    Bar colour  = peptide count at each position (colour-coded via cmap).

    Uses pre-processed data populated by process_available_data() so all
    filter-type / data-source selections are already reflected.

    Returns a plotly.graph_objects.Figure, or None on failure.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        from matplotlib.colors import to_rgba
    except ImportError:
        logger.error("plotly is required for compact heatmap: pip install plotly")
        return None

    var_keys = list(available_data_variables_dict.keys())
    num_vars = len(var_keys)
    if num_vars == 0:
        return None

    max_seq_len = max(
        len(available_data_variables_dict[v].get('protein_sequence', ''))
        for v in var_keys
    )

    # ── y-axis tick values (uniform across all subplots) ──────────────────────
    tick1 = float(y_ticks[0]) if y_ticks else 0.0
    tick2 = float(y_ticks[1]) if len(y_ticks) > 1 else (tick1 + (float(y_ticks[-1]) if y_ticks else 1.0)) / 2
    tick3 = float(y_ticks[-1]) if y_ticks else 1.0

    # Pre-compute shared y-axis scale (identical for every subplot row). A log
    # transform transforms the plotted bar HEIGHTS (below) and the axis is then
    # plain linear over log10 / log2 values — the ticks read as log values, not
    # raw abundance, so log10 vs log2 are visibly distinct.
    _do_log = bool(log_transform) and tick1 > 0
    _lsub = _log_sub(log_base)
    if _do_log:
        _lb_fn = _log_base_fn(log_base)
        _lt1, _lt2, _lt3 = (float(_lb_fn(max(tick1, 1e-300))),
                            float(_lb_fn(max(tick2, 1e-300))),
                            float(_lb_fn(max(tick3, 1e-300))))
        _y_range   = [_lt1, _lt3]
        _axis_type = 'linear'
        _tick_vals = [_lt1, _lt2, _lt3]
        _tick_text = [f"{v:.2f}" for v in _tick_vals]
    else:
        _y_range   = [tick1, tick3]
        _axis_type = 'linear'
        _tick_vals = [tick1, tick2, tick3]
        _tick_text = [f"{v:.3g}" for v in _tick_vals]

    # ── build subplots ─────────────────────────────────────────────────────────
    # Plotly caps vertical_spacing at 1/(rows-1), so a fixed 0.03 fails once
    # enough proteins are selected to push past ~34 rows. Shrink to fit.
    _v_spacing = min(0.03, 0.8 / (num_vars - 1)) if num_vars > 1 else 0
    fig = make_subplots(
        rows=num_vars,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=_v_spacing,
        subplot_titles=[""] * num_vars,
    )

    def _rgba_str(mpl_color):
        """Convert matplotlib colour (rgba tuple or 'white') to Plotly rgba string."""
        if mpl_color == 'white':
            return 'rgba(255,255,255,1.0)'
        try:
            r, g, b, a = to_rgba(mpl_color)
            return f'rgba({int(r * 255)},{int(g * 255)},{int(b * 255)},{a:.3f})'
        except Exception:
            return 'rgba(200,200,200,1.0)'

    # ── per-variable bar traces ────────────────────────────────────────────────
    for i, var_key in enumerate(var_keys):
        var_data = available_data_variables_dict[var_key]
        row = i + 1
        label = var_data.get('label', var_key)

        # prefer pre-processed data written by process_available_data()
        counts_raw = var_data.get('peptide_counts')
        ms_data    = var_data.get('ms_data_list', [])
        aa_list    = var_data.get('AA_list', [])

        if counts_raw is None or len(ms_data) == 0:
            # fall back to raw heatmap_df
            heatmap_df = var_data.get('heatmap_df')
            if heatmap_df is None or heatmap_df.empty:
                continue
            counts_list = heatmap_df['count'].tolist()
            ms_data     = heatmap_df['average'].fillna(0).tolist()
            aa_list     = heatmap_df['AA'].tolist()
        else:
            counts_list = list(counts_raw)

        positions = list(range(1, len(counts_list) + 1))

        # colour bars by peptide count
        bar_colors_mpl = get_grouped_colors(counts_list, max_count, num_colors, plot_zero, cmap)
        bar_colors = [_rgba_str(c) for c in bar_colors_mpl]

        # bar heights — keep a tiny floor so log scale stays happy
        min_height = max(tick1 * 0.01, 1e-12) if tick1 > 0 else 1e-10
        raw_heights = [
            max(float(ms), min_height) if (ms and float(ms) > 0) else min_height
            for ms in ms_data
        ]
        if _do_log:
            # Transform the plotted heights; sub-floor values sit at the bottom
            # of the shared range so they read as ~zero, same as the linear view.
            heights = [
                (_log_value(rh, log_base) if rh > min_height else _y_range[0])
                for rh in raw_heights
            ]
            heights = [h if h is not None else _y_range[0] for h in heights]
        else:
            heights = raw_heights

        fig.add_trace(
            go.Bar(
                x=positions,
                y=heights,
                marker_color=bar_colors,
                marker_line_width=0,
                width=1.0,
                hoverinfo='text',
                hovertext=[
                    f"Position: {pos}<br>Amino Acid: {aa}<br>"
                    f"Peptide Count: {cnt}<br>Averaged Abundance: {rh:.3e}"
                    + (f"<br>Log<sub>{_lsub}</sub> Abundance: {h:.2f}" if _do_log else "")
                    for pos, aa, cnt, rh, h in zip(positions, aa_list, counts_list, raw_heights, heights)
                ],
                showlegend=False,
                name=label,
            ),
            row=row, col=1,
        )

        # y-axis settings — uniform scale enforced via autorange=False; labels in legend only
        fig.update_yaxes(
            title=dict(text=label, font=dict(size=font_size_yaxis_label), standoff=5),
            type=_axis_type,
            range=_y_range,
            autorange=False,
            tickvals=_tick_vals,
            showticklabels=False,
            ticks='',
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(0,0,0,0.3)',
            # Uniform scale across subplots is the whole point of this view, and
            # the tick values are published in the legend — so zoom x only.
            fixedrange=True,
            row=row, col=1,
        )

    # ── enforce uniform y-axis scale on every subplot (global pass) ───────────
    fig.update_yaxes(
        type=_axis_type,
        range=_y_range,
        autorange=False,
        fixedrange=True,
        tickvals=_tick_vals,
        showticklabels=False,
        ticks='',
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(0,0,0,0.3)',
    )

    # ── legend: peptide count colour scale ────────────────────────────────────
    # Headers are native Plotly legend-group titles (not bold placeholder items),
    # so they align flush-left and share one grouptitlefont — matching the
    # landscape renderer and the Data Analysis legend-title convention.
    legend_title_pep  = legend_title[1] if len(legend_title) > 1 else 'Peptide Counts:'
    legend_title_abun = legend_title[2] if len(legend_title) > 2 else 'Average Abundance:'
    if _do_log:
        legend_title_abun = f'Log<sub>{_lsub}</sub> {legend_title_abun}'

    # Legend placement — 'right' keeps one combined legend beside the plot,
    # 'below' gives each group its own unit under the x-axis, one per row.
    # See LEGEND_POSITIONS.
    legend_position = normalize_legend_position(legend_position)
    legend_below = legend_position == 'below'
    _legend_kwargs, _legend_units = _make_legend_placer(legend_position, font_size_legend)

    heatmap_legend_handles, heatmap_legend_labels = create_heatmap_legend_handles(
        cmap, num_colors, max_count, plot_zero
    )
    for lbl, handle in zip(heatmap_legend_labels, heatmap_legend_handles):
        fc = handle.get_facecolor()
        r, g, b, a = fc
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(
                size=14, color=f'rgba({int(r*255)},{int(g*255)},{int(b*255)},{a:.3f})',
                symbol='square', line=dict(width=1, color='black'),
            ),
            name=lbl,
            showlegend=True,
            **_legend_kwargs(legend_title_pep, 'pep_count', entry=lbl),
        ))

    # ── legend: abundance y-range reference ───────────────────────────────────
    for tick_disp, tick_lbl in zip(_tick_text, ['Min', 'Mid', 'Max']):
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(size=10, symbol='line-ew', color='black', line=dict(width=2)),
            name=f"{tick_lbl}: {tick_disp}",
            showlegend=True,
            **_legend_kwargs(legend_title_abun, 'abundance_ticks',
                                     entry=f"{tick_lbl}: {tick_disp}"),
        ))

    # ── x-axes ────────────────────────────────────────────────────────────────
    x_label_str = xaxis_label or 'Protein Sequence Position'
    for i in range(1, num_vars + 1):
        common = dict(
            showgrid=True,
            gridwidth=1,
            gridcolor='rgba(200,200,200,0.5)',
            range=[-1, max_seq_len + 0.5],
            autorange=False,
            minallowed=-1,
            maxallowed=max_seq_len + 0.5,
            zeroline=False,
        )
        if i == num_vars:
            fig.update_xaxes(
                title=dict(text=x_label_str, font=dict(size=font_size_xaxis_label)),
                dtick=20, tickfont=dict(size=font_size_xaxis_tick),
                **common,
                row=i, col=1,
            )
        else:
            fig.update_xaxes(
                showticklabels=False,
                **common,
                row=i, col=1,
            )

    # ── overall layout ────────────────────────────────────────────────────────
    # No default title (see landscape renderer): the "<name> Sequence" x-axis and
    # the page heading already name the protein. Optional user title only.
    title_str = (plot_title or '').strip()

    row_height = max(130, 600 // max(num_vars, 1))

    # ── margins + legend geometry (see the landscape renderer for the rationale) ──
    top_margin   = 60 if title_str else 30
    right_margin = 40 if legend_below else 180
    base_height  = max(row_height * num_vars + 80, 300)
    legend_meta  = None

    if legend_below:
        legend_layout, fig_height, bottom_margin, legend_meta = _below_legend_layout(
            _below_legend_plan(_legend_units, font_size_legend),
            True, top_margin, base_height, font_size_legend)
    else:
        fig_height    = base_height
        bottom_margin = LEGEND_BASE_BOTTOM_PX
        legend_layout = dict(legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.01,
            bordercolor="black",
            borderwidth=1,
            orientation="v",
            bgcolor="rgba(255,255,255,0.9)",
            # One shared font size for both legend groups' headers AND items —
            # Peptide Counts and Average Abundance both live in this one legend.
            font=dict(size=font_size_legend),
            grouptitlefont=dict(size=font_size_legend, color='black'),
            itemsizing='constant',
            groupclick="toggleitem",
        ))

    if legend_meta is not None:
        # See the landscape renderer: lets the browser / export re-pack the
        # units once the real plotting width is known.
        fig.update_layout(meta=legend_meta)

    fig.update_layout(
        title=dict(text=title_str, font=dict(size=font_size_plot_title, color='black'), x=0.5),
        font=dict(family=PLOTLY_FONT_FAMILY, color='black', size=12),
        height=fig_height,
        bargap=0,
        bargroupgap=0,
        plot_bgcolor='white',
        paper_bgcolor='white',
        showlegend=True,
        margin=dict(l=60, r=right_margin, t=top_margin, b=bottom_margin),
        **legend_layout,
    )

    return fig


"""_______________________________Portrait Plot______________________________________________"""
