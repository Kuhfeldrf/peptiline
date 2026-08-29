"""
Plotting service for the Data Analysis web app.
Extracted from data_analysis.ipynb Plotter class.
All functions take a DataAnalysisState object and return a Plotly Figure.
"""
import traceback

import numpy as np
import pandas as pd

from .data_processor import (
    DataAnalysisState,
    contains_function,
    get_color_sequence,
    get_single_color,
    redact_string_descriptions,
)

# Largest square edge, in px, for the correlation scatter-plot matrix. Keeps a
# many-group SPLOM inside the browser viewport instead of overflowing the page.
MAX_SPLOM_SIZE = 1000

# ── Shared figure styling ──────────────────────────────────────────────────
# Single source of truth so every Data Analysis plot type (bar, stacked, pie,
# correlation, SPLOM) shares one typeface, hover style and margins. Plotly's
# layout.font.family cascades to titles, axes, legend and annotations — which
# set only size/colour, not family — so defining the global font once here
# standardises the type across the whole figure. Change these to restyle all
# figures at once.
FONT_FAMILY = 'Arial, Helvetica, sans-serif'
BASE_FONT   = dict(family=FONT_FAMILY, color='black')
HOVERLABEL  = dict(bgcolor='white', font_size=12, font_family=FONT_FAMILY)
PLOT_MARGIN = dict(t=100, l=100, r=100, b=100)

# Below this max plotted value, an abundance chart defaults to a linear y-axis
# instead of the usual log scale (the log-transform checkbox still overrides).
LINEAR_SCALE_THRESHOLD = 1e4

# Shown when the user enables Show Stats on a Relative (composition) plot. The
# group comparison is defined on the replicate-level absolute values; relative
# percentages are compositional (constrained to sum to 100%), so their replicate
# variance is not the replicate SEM and an ANOVA on them would violate the
# independence assumption. Switch to Absolute (or Peptide Count) for significance.
RELATIVE_STATS_NOTE = (
    'Significance not shown for relative (composition) metrics — the group '
    'comparison is defined on absolute values. Switch to an Absolute metric.'
)


def _log_arr(v, base=10):
    """Base-aware elementwise log (base 10 or 2) for array-likes / scalars."""
    return np.log2(v) if int(base) == 2 else np.log10(v)


def _log_sub(base):
    """Subscript string for axis labels: 'Log<sub>10</sub>' / 'Log<sub>2</sub>'."""
    return '2' if int(base) == 2 else '10'


def _log_hover_line(state, abs_value):
    """Hover fragment showing the log-transformed value when a log transform is
    active — so the tooltip reports the value as plotted, not just raw abundance.
    Returns '' (no extra line) when log transform is off or the value is
    non-positive."""
    if not getattr(state, 'log_transform', False):
        return ''
    lv = _safe_log(abs_value, getattr(state, 'log_base', 10))
    if lv is None:
        return ''
    return f"Log<sub>{_log_sub(state.log_base)}</sub> {state.metric_name}: {lv:.2f}<br>"


def _safe_log(val, base=10):
    try:
        v = float(val)
    except (TypeError, ValueError):
        return None
    if np.isnan(v) or v <= 0:
        return None
    return float(_log_arr(v, base))


def _display_sem(sem, value, use_log, base=10):
    """Transform a raw replicate SEM so it matches the plotted y value.

    Linear axis: return the SEM unchanged. Log axis: apply the delta-method
    transform SEM/(value·ln(base)) (the same conversion the totals plot uses at
    ``plot_total_peptides``). Returns 0.0 for a non-positive SEM or, under log,
    a non-positive value — which renders as no error bar.
    """
    try:
        s = float(sem)
    except (TypeError, ValueError):
        return 0.0
    if s <= 0 or np.isnan(s):
        return 0.0
    if use_log:
        try:
            v = float(value)
        except (TypeError, ValueError):
            return 0.0
        return s / (v * np.log(int(base))) if v > 0 else 0.0
    return s


# ---------------------------------------------------------------------------
# Title / label helpers
# ---------------------------------------------------------------------------

def _make_title(state: DataAnalysisState, kind: str = 'bar') -> str:
    """Build the dynamic plot title from the active selections.

    Bar plots: "<Absolute|Relative> <Abundance|Count> Distribution[ by <Protein|
    Function>]". The metric qualifier states absolute vs relative; the orientation
    qualifier is kept only when it is meaningful (By Protein / By Function) and
    dropped for By Sample.

    Correlation plots use a distinct, correlation-specific title (they show a
    pairwise relationship, not a distribution). A user-supplied ``plot_title``
    always overrides.
    """
    if state.plot_title:
        return state.plot_title

    if kind == 'correlation':
        # Correlation is computed on absolute abundance (Avg_ columns); relative
        # composition isn't meaningful for a pairwise correlation, so no metric
        # qualifier is needed here.
        return 'Pairwise Abundance Correlation'

    metric_qual = 'Relative' if state.is_relative else 'Absolute'
    base_metric = 'Count' if state.use_count else 'Abundance'
    title = f'{metric_qual} {base_metric} Distribution'
    if state.orientation in ('By Protein', 'By Function'):
        title += ' ' + state.orientation.replace('By ', 'by ')
    return title


# ---------------------------------------------------------------------------
# Plot 1: Total peptides (By Sample, No Filter / Both)
# ---------------------------------------------------------------------------

def _plot_totals_relative(state: DataAnalysisState):
    """Relative version of the sample-totals chart (No Filter / Both, By Sample):
    each sample as a percentage of the grand total across the selected samples.

    Relative composition is descriptive only — no SEM error bars, no log axis, and
    no significance markers (those are defined on absolute values; see
    RELATIVE_STATS_NOTE). This is the counterpart that was missing, so toggling
    Relative under No Filter / Both had no effect and the plot stayed absolute.
    """
    import plotly.graph_objects as go

    data = state.total_peptide_results_dict
    if not data:
        return None

    groups = list(data.keys())
    key = 'unique_peptides' if state.use_count else 'total_Abundance'
    total = sum(data[g][key] for g in groups) or 1.0
    rels = [data[g][key] / total * 100.0 for g in groups]
    first_color = get_single_color(state.color_scheme)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=groups, y=rels, name='Relative ' + state.metric_name,
        marker=dict(color=first_color, line=dict(color='black', width=1)),
        hovertext=[
            f"Group: {g}<br>Contribution: {r:.2f}%<br>"
            f"{state.metric_name}: {data[g][key]:{state.num_format}}"
            for g, r in zip(groups, rels)
        ],
        hoverinfo='text',
    ))
    fig.add_trace(go.Scatter(
        x=groups, y=[r + 1.5 for r in rels], mode='text',
        text=[f"{r:.1f}%" for r in rels], textposition='top center',
        textfont=dict(size=state.font_size_value_label, color='black'),
        showlegend=False, hoverinfo='none',
    ))

    y_title = 'Relative ' + state.metric_name
    yaxis_kw = _axis_style(state.font_size_ylabel, state.font_size_ytick)
    yaxis_kw['range'] = [0, 100]
    yaxis_kw['tickformat'] = '.1f'
    fig.update_layout(**{
        **_common_layout(state),
        'showlegend': False,
        'hoverlabel': HOVERLABEL,
        'xaxis': dict(title=state.xlabel or '', tickangle=45, **_axis_style(state.font_size_xlabel, state.font_size_xtick)),
        'yaxis': dict(title=state.ylabel or y_title, **yaxis_kw),
    })
    return fig


def plot_total_peptides(state: DataAnalysisState):
    import plotly.graph_objects as go

    data = state.total_peptide_results_dict
    if not data:
        return None

    # Relative metric: each sample as a share of the grand total. Without this the
    # Relative toggle was ignored for No Filter / Both (this chart only rendered
    # absolute values). See _plot_totals_relative.
    if state.is_relative:
        return _plot_totals_relative(state)

    first_color = get_single_color(state.color_scheme)
    use_log = state.log_transform
    groups = list(data.keys())

    abundances = [data[g]['total_Abundance'] for g in groups]
    counts = [data[g]['unique_peptides'] for g in groups]
    abundance_sems = [data[g]['abundance_sem'] for g in groups]
    count_sems = [data[g]['count_sem'] for g in groups]

    if use_log:
        base = state.log_base
        abundances = [_safe_log(v, base) for v in abundances]
        counts = [_safe_log(v, base) for v in counts]
        skipped_log_zero = sum(1 for v in (abundances if not state.use_count else counts) if v is None)
        if skipped_log_zero:
            state.warnings.append(
                f'{skipped_log_zero} group(s) had a zero or missing total and were left blank '
                'under log transform (log of zero/negative is undefined). Turn off log '
                'transform to display them.'
            )
        abundance_sems = [
            data[g]['abundance_sem'] / (data[g]['total_Abundance'] * np.log(base))
            if data[g]['total_Abundance'] > 0 else 0
            for g in groups
        ]
        count_sems = [
            data[g]['count_sem'] / (data[g]['unique_peptides'] * np.log(base))
            if data[g]['unique_peptides'] > 0 else 0
            for g in groups
        ]
        y_prefix = f'Log<sub>{_log_sub(base)}</sub> '
        tickfmt = '.2f'
    else:
        y_prefix = ''
        tickfmt = '.1e'

    if state.use_count:
        count_label_text = [('' if c is None else f"{c:.2f}") if use_log else f"{int(c):,}" for c in counts]
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=groups, y=counts, name='Peptide Count',
            marker=dict(color=first_color, line=dict(color='black', width=1)),
            error_y=dict(type='data', array=count_sems, visible=True,
                         thickness=1.5, width=4, color='#000000'),
            hovertemplate=(
                "Group: %{x}<br>"
                f"{y_prefix}Unique Peptides: %{{y:{'.2f' if use_log else '.0f'}}}<br>"
                f"SEM: %{{error_y.array:{'.2f' if use_log else '.1f'}}}<br>"
                "<extra></extra>"
            ),
        ))
        fig.add_trace(go.Scatter(
            x=groups,
            y=[(0 if c is None else c) + (s * 1.2) for c, s in zip(counts, count_sems)],
            mode='text',
            text=count_label_text,
            textposition='top center',
            textfont=dict(size=state.font_size_value_label),
            showlegend=False,
            hoverinfo='none',
        ))
        y_axis_title = f'{y_prefix}Unique Peptide Count'
        if use_log:
            yaxis_kw = dict(tickformat=tickfmt, **_axis_style(state.font_size_ylabel, state.font_size_ytick))
        else:
            yaxis_kw = dict(tickformat=',d', exponentformat='none',
                            **_axis_style(state.font_size_ylabel, state.font_size_ytick))
        fig.update_layout(
            **_common_layout(state),
            showlegend=False,
            hoverlabel=HOVERLABEL,
            # Totals plot is always by sample and the group ticks are self-
            # explanatory, so no default axis title — only what the user supplies.
            xaxis=dict(title=state.xlabel or '', tickangle=45, **_axis_style(state.font_size_xlabel, state.font_size_xtick)),
            yaxis=dict(title=y_axis_title, **yaxis_kw),
        )
    else:
        abundance_label_text = [('' if a is None else f"{a:.2f}") if use_log else f"{a:.2e}" for a in abundances]
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=groups, y=abundances, name='Abundance',
            marker=dict(color=first_color, line=dict(color='black', width=1)),
            error_y=dict(type='data', array=abundance_sems, visible=True,
                         thickness=1.5, width=4, color='#000000'),
            hovertemplate=(
                "Group: %{x}<br>"
                f"{y_prefix}Total Abundance: %{{y:{'.2f' if use_log else '.2e'}}}<br>"
                f"SEM: %{{error_y.array:{'.2f' if use_log else '.2e'}}}<br>"
                "<extra></extra>"
            ),
        ))
        fig.add_trace(go.Scatter(
            x=groups,
            y=[(0 if a is None else a) + s for a, s in zip(abundances, abundance_sems)],
            mode='text',
            text=abundance_label_text,
            textposition='top center',
            textfont=dict(size=state.font_size_value_label),
            showlegend=False,
            hoverinfo='none',
        ))
        y_axis_title = f'{y_prefix}Summed Abundance'
        if use_log:
            yaxis_kw = dict(tickformat=tickfmt, **_axis_style(state.font_size_ylabel, state.font_size_ytick))
        elif max(abundances, default=0) < LINEAR_SCALE_THRESHOLD:
            # Small values read better on a linear axis than a log axis; the
            # log-transform checkbox is untouched and still forces log when ticked.
            # Below this threshold (including abundance data that was already
            # log-transformed upstream) plain integers read better than scientific notation.
            yaxis_kw = dict(type='linear', tickformat=',d', exponentformat='none',
                            **_axis_style(state.font_size_ylabel, state.font_size_ytick))
        else:
            yaxis_kw = dict(type='log', exponentformat='e',
                            **_axis_style(state.font_size_ylabel, state.font_size_ytick))
        fig.update_layout(
            **_common_layout(state),
            showlegend=False,
            hoverlabel=HOVERLABEL,
            xaxis=dict(title=state.xlabel or '', tickangle=45, **_axis_style(state.font_size_xlabel, state.font_size_xtick)),
            yaxis=dict(title=state.ylabel or y_axis_title, **yaxis_kw),
        )

    if state.show_significance and not state.is_relative:
        if state.use_count:
            tops = {g: (0 if counts[i] is None else counts[i]) + count_sems[i] for i, g in enumerate(groups)}
        else:
            tops = {g: (0 if abundances[i] is None else abundances[i]) + abundance_sems[i] for i, g in enumerate(groups)}
        # Log-scale axis is used only for abundance, without the log checkbox, and
        # above the linear-scale threshold (see the yaxis block above) — the marker
        # placement must match whichever axis actually rendered.
        is_log_axis = (not state.use_count) and (not use_log) \
            and max(abundances, default=0) >= LINEAR_SCALE_THRESHOLD
        _add_totals_significance(fig, groups, state.total_reps_dict,
                                 'count' if state.use_count else 'abundance', tops, is_log_axis,
                                 method=state.significance_method,
                                 sig_font_size=state.font_size_value_label)
    return fig


# ---------------------------------------------------------------------------
# Plot 2: Grouped bar plot
# ---------------------------------------------------------------------------

def _significance_caption(method_label, excluded, min_n, tested_any, reason=None):
    """Build the footnote that keeps the significance overlay honest.

    When a comparison ran, it names the test + alpha, lists any groups dropped for
    too few replicates, and — because power is dire at these replicate counts —
    warns that a non-significant result is inconclusive when n is small. When
    nothing could be tested, it says *why* rather than leaving the user staring at
    a checkbox that silently did nothing.
    """
    if not tested_any:
        why = reason or 'need at least 3 replicates in at least 2 groups'
        return f'Significance ({method_label}) not shown: {why}.'
    parts = [f'{method_label}, α = 0.05']
    if excluded:
        parts.append('excluded (fewer than 3 replicates): ' + ', '.join(excluded))
    if min_n is not None and min_n <= 4:
        parts.append(
            f'smallest group n = {min_n}: power is low, so a shared letter / "ns" '
            'means not resolvable at this sample size, not "no difference"')
    return '; '.join(parts) + '.'


def _add_significance_caption(fig, text):
    """Append the significance footnote below the x-axis.

    Rendered as an extra line of the x-axis title with ``automargin=True`` (rather
    than a fixed paper-y annotation) so Plotly places it *below* the tick labels
    no matter how long or steeply-rotated the category labels are — a fixed-offset
    annotation collides with long rotated labels on the By-Function / By-Protein
    charts. The footnote is kept visually distinct (smaller, grey) from any real
    axis title it sits under.
    """
    if not text:
        return
    try:
        existing = fig.layout.xaxis.title.text or ''
    except Exception:
        existing = ''
    caption = f'<span style="font-size:11px;color:#555555">{text}</span>'
    # Strip a caption already appended on a prior pass (defensive; each figure is
    # built fresh, so this normally just guards against double-calls).
    existing = existing.split('<br> <br>', 1)[0]
    new_title = f'{existing}<br> <br>{caption}' if existing.strip() else caption
    fig.update_xaxes(title_text=new_title, automargin=True)


def _offset_fn(tops_vals, is_log_axis):
    """Build ``above(base, steps)`` placing markers ``steps`` gaps above a value.

    The gap is sized to the *visible* spread of the bars so markers clear labels
    without ballooning: on a linear axis it's a fraction of the value span; on a
    log axis it's a fraction of the log-range of the bar tops (a fixed
    multiplicative factor would swallow a sub-decade view). One gap ≈ clears a
    line of label text.
    """
    import math

    vals = [v for v in tops_vals if v is not None]
    if is_log_axis:
        logs = [math.log10(v) for v in vals if v > 0] or [0.0]
        step = 0.11 * max(max(logs) - min(logs), 0.3)  # decades per gap

        def above(base, steps):
            return 10 ** (math.log10(base) + steps * step) if base and base > 0 else base
    else:
        span = max(vals, default=0) or 1.0
        step = 0.06 * span

        def above(base, steps):
            return (base or 0) + steps * step
    return above


def _add_totals_significance(fig, groups, total_reps, metric_key, tops, is_log_axis,
                             method='tukey', x_categorical=True, sig_font_size=None):
    """Significance markers for a chart whose x-clusters are the sample groups.

    Used both by the totals bar chart and by grouped bars in "By Sample"
    orientation (where each sample is a cluster). Compares the sample groups on
    their replicate totals: two groups get a bracket, three or more a compact-letter
    display.

    ``tops`` maps each group to the data-space top of its bar+SEM. Markers are
    lifted a gap (``above(base, 1)``) to clear the numeric value label the totals
    chart prints, and the two-group bracket is drawn in data coordinates at the max
    of the two bars — never a fixed height that could collide with an error bar or
    label. Markers are emitted as Scatter traces so the y-axis auto-extends on both
    linear and log axes.

    ``x_categorical`` selects how the scatter-based markers reference x: True for a
    true Plotly category axis (totals plot — reference by group NAME; a numeric x
    would spawn stray "0"/"1" categories), False for a numeric axis with
    tickvals/ticktext (grouped By-Sample — reference by position index). Shapes and
    annotations always use the numeric position index, which is valid on both.
    ``method`` selects the test; see services/stats.py.
    """
    import plotly.graph_objects as go

    from . import stats

    rep_arrays = {g: list(total_reps.get(g, {}).get(metric_key, []) or []) for g in groups}
    res = stats.compare_groups(rep_arrays, method=method)
    excluded = [redact_string_descriptions(g) for g in res.get('excluded', [])]
    if not res['testable']:
        _add_significance_caption(fig, _significance_caption(
            res['method'], excluded, None, False, res.get('reason')))
        return

    above = _offset_fn(tops.values(), is_log_axis)
    idx = {g: i for i, g in enumerate(groups)}
    # One gap above each bar clears its numeric value label; markers go above that.
    label_top = {g: above(tops.get(g, 0), 1) for g in groups}
    anchor_y = None
    star_size = sig_font_size if sig_font_size is not None else 14
    letter_size = sig_font_size if sig_font_size is not None else 13

    if len(groups) == 2:
        g1, g2 = groups
        stars = stats.significance_stars(res['pairwise'].get(frozenset((g1, g2))))
        y = above(max(label_top[g1], label_top[g2]), 1)
        # Bracket via shapes/annotation: on a category axis, numeric x in a SHAPE
        # is a position index (0, 1, 0.5 = midpoint) — unlike a scatter trace,
        # where numbers would create spurious categories. Each riser starts above
        # its own bar's value label; the crossbar sits above the taller of the two.
        line = dict(color='black', width=1)
        fig.add_shape(type='line', xref='x', yref='y', x0=idx[g1], x1=idx[g1],
                      y0=label_top[g1], y1=y, line=line)
        fig.add_shape(type='line', xref='x', yref='y', x0=idx[g2], x1=idx[g2],
                      y0=label_top[g2], y1=y, line=line)
        fig.add_shape(type='line', xref='x', yref='y', x0=idx[g1], x1=idx[g2],
                      y0=y, y1=y, line=line)
        # Plotly quirk: layout ANNOTATIONS on a log axis take log10(value) for y
        # (unlike shapes and scatter traces, which use the raw value). Without this
        # the star lands at 10**y and is clipped off the top of the chart.
        star_y = float(np.log10(y)) if (is_log_axis and y > 0) else y
        fig.add_annotation(xref='x', yref='y', x=(idx[g1] + idx[g2]) / 2, y=star_y,
                           text=stars, showarrow=False, yanchor='bottom',
                           font=dict(size=star_size, color='black'))
        anchor_y = above(y, 1.5)
    else:
        # Compact-letter display as a scatter-text trace, auto-extends the y-range.
        letter_groups = [g for g in groups if res['letters'].get(g)]
        fig.add_trace(go.Scatter(
            x=(letter_groups if x_categorical else [idx[g] for g in letter_groups]),
            y=[label_top[g] for g in letter_groups],
            mode='text',
            text=[f"<b>{res['letters'][g]}</b>" for g in letter_groups],
            textposition='top center',
            textfont=dict(size=letter_size, color='black'),
            showlegend=False, hoverinfo='none',
        ))
        anchor_y = above(max((label_top[g] for g in letter_groups),
                             default=max(tops.values(), default=1.0)), 1.5)

    # Shapes/annotations do not auto-extend the axis, so drop an invisible marker to
    # pull the y-range up and keep the markers on-canvas (linear and log axes).
    if anchor_y is not None:
        fig.add_trace(go.Scatter(
            x=([groups[0]] if x_categorical else [0]), y=[anchor_y], mode='markers',
            marker=dict(opacity=0), showlegend=False, hoverinfo='none',
        ))

    _add_significance_caption(fig, _significance_caption(
        res['method'], excluded, res.get('min_n'), True))


def _add_grouped_significance(fig, categories, group_order, geom, reps_lookup, metric_key,
                              method='tukey', sig_font_size=None):
    """Overlay ANOVA + Tukey HSD significance markers on a grouped bar chart.

    For each category cluster (one function or one protein) the sample groups are
    compared on their replicate-level values. Two groups get a significance
    bracket (``*``/``**``/``***``/``ns``); three or more get a compact-letter
    display (bars sharing a letter are not significantly different). Markers are
    only drawn where a valid test exists (≥3 replicates per group); the y-range
    is extended so they sit clear of the tallest bar.

    Parameters
    ----------
    geom : dict          {(category_index, group): {'x': centre, 'top': bar+SEM top}}
    reps_lookup : callable   category name -> {group: {'abundance': [...], 'count': [...]}}
    metric_key : str     'abundance' or 'count'
    """
    from . import stats

    tops = [v['top'] for v in geom.values() if v and v['top'] is not None]
    axis_max = max(tops, default=0.0)
    off = (axis_max * 0.04) if axis_max > 0 else 1.0
    needed_top = axis_max
    any_marker = False
    method_label = stats._METHOD_LABELS.get(
        stats._METHOD_ALIASES.get(str(method).strip().lower(), 'tukey'))
    excluded, min_ns, last_reason = set(), [], None
    star_size = sig_font_size if sig_font_size is not None else 13
    letter_size = sig_font_size if sig_font_size is not None else 12

    for ci, cat in enumerate(categories):
        reps = reps_lookup(cat) or {}
        rep_arrays = {g: list(reps.get(g, {}).get(metric_key, []) or []) for g in group_order}
        res = stats.compare_groups(rep_arrays, method=method)
        excluded.update(res.get('excluded', []))
        if not res['testable']:
            last_reason = last_reason or res.get('reason')
            continue
        if res.get('min_n') is not None:
            min_ns.append(res['min_n'])
        cat_geom = {g: geom.get((ci, g)) for g in group_order if geom.get((ci, g))}
        if not cat_geom:
            continue
        any_marker = True
        cluster_top = max(v['top'] for v in cat_geom.values())

        # Two displayed groups → significance bracket with stars. Each vertical
        # riser starts just above its OWN bar+SEM top and the horizontal sits above
        # the taller of the two, so the bracket never overlaps either error bar.
        if len(group_order) == 2 and all(g in cat_geom for g in group_order):
            g1, g2 = group_order
            label = stats.significance_stars(res['pairwise'].get(frozenset((g1, g2))))
            t1, t2 = cat_geom[g1]['top'], cat_geom[g2]['top']
            x1, x2 = cat_geom[g1]['x'], cat_geom[g2]['x']
            y = max(t1, t2) + off * 1.5
            fig.add_shape(type='line', x0=x1, x1=x1, y0=t1 + off * 0.4, y1=y,
                          line=dict(color='black', width=1))
            fig.add_shape(type='line', x0=x2, x1=x2, y0=t2 + off * 0.4, y1=y,
                          line=dict(color='black', width=1))
            fig.add_shape(type='line', x0=x1, x1=x2, y0=y, y1=y, line=dict(color='black', width=1))
            fig.add_annotation(x=(x1 + x2) / 2, y=y + off * 0.2, text=label, showarrow=False,
                               font=dict(size=star_size, color='black'), yanchor='bottom')
            needed_top = max(needed_top, y + off * 2)
        else:
            # Three or more groups → compact-letter display above each bar.
            for g, letter in res['letters'].items():
                cg = cat_geom.get(g)
                if cg and letter:
                    fig.add_annotation(x=cg['x'], y=cg['top'] + off * 0.35, text=f"<b>{letter}</b>",
                                       showarrow=False, font=dict(size=letter_size, color='black'),
                                       yanchor='bottom')
            needed_top = max(needed_top, cluster_top + off * 2)

    if any_marker and needed_top > axis_max:
        # Linear grouped-bar axis (log-transform pre-logs the values), so a plain
        # [0, top] range is valid and leaves headroom for the markers.
        fig.update_yaxes(range=[0, needed_top + off])

    excluded_disp = [redact_string_descriptions(g) for g in group_order if g in excluded]
    _add_significance_caption(fig, _significance_caption(
        method_label, excluded_disp, min(min_ns) if min_ns else None,
        any_marker, last_reason))


def create_grouped_bar_plot(state: DataAnalysisState):
    import plotly.graph_objects as go

    orientation = state.orientation
    plot_filter = state.plot_filter
    selected_groups = state.selected_groups
    use_count = state.use_count
    use_log = state.log_transform
    log_base = state.log_base
    metric_key = 'count' if use_count else 'abundance'
    show_sig = state.show_significance and not state.is_relative

    if orientation == 'By Function' or (orientation == 'By Sample' and plot_filter in ('Selected Function(s)', 'Functional vs Non-Functional Peptides')):
        if not state.function_distribution_dict:
            return None
        if orientation == 'By Sample':
            categories = selected_groups
            bar_groups = state.selected_functions
            color_seq = get_color_sequence(len(bar_groups), state.color_scheme)
            color_map = {item: color_seq[i] for i, item in enumerate(bar_groups)}
        else:
            categories = state.selected_functions
            bar_groups = selected_groups
            color_seq = get_color_sequence(len(selected_groups), state.color_scheme)
            color_map = {g: color_seq[i] for i, g in enumerate(selected_groups)}

        # Other Functions always shown in grey
        if 'Other Functions' in color_map:
            color_map['Other Functions'] = '#808080'

        display_cats = [redact_string_descriptions(c) for c in categories]
        n_bars = len(bar_groups)
        bar_width = 0.8 / max(n_bars, 1)
        fig = go.Figure()
        sig_geom = {}  # {(category_index, group): {'x': centre, 'top': bar+SEM top}}
        bysample_tops = {}  # By Sample: {sample: max bar+SEM in its cluster}
        all_values = []
        skipped_log_zero = 0

        for idx, bar_group in enumerate(bar_groups):
            x_pos = [i + (idx - n_bars / 2 + 0.5) * bar_width for i in range(len(categories))]
            values, hover, sems = [], [], []
            disp_bg = redact_string_descriptions(bar_group)

            for ci_cat, cat in enumerate(categories):
                if orientation == 'By Sample':
                    item, group = bar_group, cat
                else:
                    item, group = cat, bar_group

                entry = state.function_distribution_dict.get(item, {})
                abs_value = entry.get('counts' if use_count else 'Abundance', {}).get(group, 0)
                rel_value = entry.get('count_relative' if use_count else 'abundance_relative', {}).get(group, 0)
                sem_raw = entry.get('count_sem' if use_count else 'abundance_sem', {}).get(group, 0)
                if state.is_relative:
                    v = rel_value
                else:
                    v = abs_value
                    if use_log:
                        v = _safe_log(v, log_base)
                        if v is None:
                            skipped_log_zero += 1
                values.append(v)
                disp_sem = _display_sem(sem_raw, abs_value, use_log, log_base)
                sems.append(disp_sem)
                # Record bar geometry for significance markers. By Function: each
                # category cluster holds one bar per sample, compared within-cluster.
                # By Sample: samples ARE the clusters, so track each cluster's top
                # and compare samples on the filtered aggregate (see below).
                if show_sig and orientation == 'By Function':
                    sig_geom[(ci_cat, group)] = {'x': x_pos[ci_cat],
                                                 'top': (v or 0) + disp_sem}
                elif show_sig and orientation == 'By Sample':
                    bysample_tops[cat] = max(bysample_tops.get(cat, 0.0), (v or 0) + disp_sem)
                item_label = 'Function' if plot_filter in ('Selected Function(s)', 'Functional vs Non-Functional Peptides') else 'Protein'
                if state.is_relative:
                    hover.append(
                        f"{item_label}: {item}<br>"
                        f"Sample: {group}<br>"
                        f"Relative Contribution: {rel_value:.1f}%<br>"
                        f"{state.metric_name}: {abs_value:{state.num_format}}"
                    )
                else:
                    hover.append(
                        f"{item_label}: {item}<br>"
                        f"Sample: {group}<br>"
                        f"{state.metric_name}: {abs_value:{state.num_format}}<br>"
                        f"{_log_hover_line(state, abs_value)}"
                        f"Relative Contribution: {rel_value:.1f}%"
                    )

            all_values.extend(values)
            fn_color = color_map.get(bar_group, '#999')
            if bar_group in ('Other Functions', 'Other Proteins'):
                fn_color = '#808080'
            # Replicate-based SEM error bars: only on absolute metrics (a percentage
            # bar's replicate SEM is not the SEM of the ratio, so it is omitted).
            bar_error_y = None
            if not state.is_relative and any(s > 0 for s in sems):
                bar_error_y = dict(type='data', array=sems, visible=True,
                                   thickness=1.2, width=3, color='#000000')
            fig.add_trace(go.Bar(
                x=x_pos, y=values, name=disp_bg,
                marker=dict(color=fn_color, line=dict(color='black', width=0.5)),
                width=bar_width * 0.9,
                error_y=bar_error_y,
                hovertext=hover, hoverinfo='text',
            ))

        if skipped_log_zero:
            state.warnings.append(
                f'{skipped_log_zero} bar(s) had a zero or missing value and were left blank '
                'under log transform (log of zero/negative is undefined). Turn off log '
                'transform to display them.'
            )
        y_title = state.ylabel or (f'Log<sub>{_log_sub(log_base)}</sub> ' if use_log else '') + state.metric_name
        finite_values = [v for v in all_values if v is not None]
        yaxis_fn = _build_grouped_yaxis(state, use_log, use_count, max(finite_values, default=0))
        x_lbl, legend_lbl = _orientation_axis_titles(state)
        fig.update_layout(
            **_common_layout(state),
            xaxis=dict(tickvals=list(range(len(categories))), ticktext=display_cats,
                       tickangle=45, title=x_lbl,
                       **_axis_style(state.font_size_xlabel, state.font_size_xtick)),
            yaxis=dict(title=y_title, **yaxis_fn),
            barmode='group',
            hoverlabel=HOVERLABEL,
            legend=_legend_style(state, legend_lbl, yanchor='top', y=1.0, xanchor='left', x=1.05),
        )
        if show_sig and orientation == 'By Function' and sig_geom:
            _add_grouped_significance(
                fig, categories, list(bar_groups), sig_geom,
                lambda fn: state.function_reps_dict.get(fn, {}), metric_key,
                method=state.significance_method, sig_font_size=state.font_size_value_label)
        elif show_sig and orientation == 'By Sample' and bysample_tops:
            # Samples are the x-clusters here, so compare the samples on their
            # replicate-level filtered totals (total_reps_dict is computed on the
            # filtered_df, i.e. the selected functions) and annotate per cluster.
            # Grouped bars use a numeric x-axis (tickvals/ticktext) -> x_categorical=False.
            _add_totals_significance(fig, categories, state.total_reps_dict, metric_key,
                                     bysample_tops, is_log_axis=False,
                                     method=state.significance_method, x_categorical=False,
                                     sig_font_size=state.font_size_value_label)
        elif state.show_significance and state.is_relative:
            _add_significance_caption(fig, RELATIVE_STATS_NOTE)
        return fig

    elif orientation == 'By Protein' or (orientation == 'By Sample' and plot_filter == 'Selected Protein(s)'):
        if state.protein_df is None or state.protein_df.empty:
            return None
        df = state.protein_df
        if orientation == 'By Sample':
            categories = selected_groups
            bar_groups_list = list(df['Description'])
            color_seq = get_color_sequence(len(bar_groups_list), state.color_scheme)
            color_map = {item: color_seq[i] for i, item in enumerate(bar_groups_list)}
        else:
            categories = list(df['Description'])
            bar_groups_list = selected_groups
            color_seq = get_color_sequence(len(selected_groups), state.color_scheme)
            color_map = {g: color_seq[i] for i, g in enumerate(selected_groups)}

        # Other Proteins always shown in grey
        if 'Other Proteins' in color_map:
            color_map['Other Proteins'] = '#808080'

        display_cats = [redact_string_descriptions(c) for c in categories]
        n_bars = len(bar_groups_list)
        bar_width = 0.8 / max(n_bars, 1)
        fig = go.Figure()
        sig_geom = {}  # {(category_index, group): {'x': centre, 'top': bar+SEM top}}
        bysample_tops = {}  # By Sample: {sample: max bar+SEM in its cluster}
        all_values = []
        skipped_log_zero = 0

        for idx, bar_group in enumerate(bar_groups_list):
            x_pos = [i + (idx - n_bars / 2 + 0.5) * bar_width for i in range(len(categories))]
            values, hover, sems = [], [], []
            disp_bg = redact_string_descriptions(bar_group)

            for ci_cat, cat in enumerate(categories):
                if orientation == 'By Sample':
                    item, group = bar_group, cat
                else:
                    item, group = cat, bar_group
                # Use Count columns when Peptide Count is selected, Avg columns for Abundance
                abs_col = f'Count_{group}' if use_count else f'Avg_{group}'
                rel_col = f'Rel_Count_{group}' if use_count else f'Rel_Avg_{group}'
                sem_col = f'SEM_Count_{group}' if use_count else f'SEM_Avg_{group}'
                row = df[df['Description'] == item]
                abs_value = float(row[abs_col].values[0]) if len(row) > 0 and abs_col in row.columns else 0
                rel_value = float(row[rel_col].values[0]) if len(row) > 0 and rel_col in row.columns else 0
                sem_raw = float(row[sem_col].values[0]) if len(row) > 0 and sem_col in row.columns else 0
                if state.is_relative:
                    v = rel_value
                else:
                    v = abs_value
                    if use_log:
                        v = _safe_log(v, log_base)
                        if v is None:
                            skipped_log_zero += 1
                values.append(v)
                disp_sem = _display_sem(sem_raw, abs_value, use_log, log_base)
                sems.append(disp_sem)
                if show_sig and orientation == 'By Protein':
                    sig_geom[(ci_cat, group)] = {'x': x_pos[ci_cat],
                                                 'top': (v or 0) + disp_sem}
                elif show_sig and orientation == 'By Sample':
                    bysample_tops[cat] = max(bysample_tops.get(cat, 0.0), (v or 0) + disp_sem)
                count_fmt = ',.0f' if use_count else state.num_format
                if state.is_relative:
                    hover.append(
                        f"Protein: {item}<br>"
                        f"Sample: {group}<br>"
                        f"Relative Contribution: {rel_value:.1f}%<br>"
                        f"{state.metric_name}: {abs_value:{count_fmt}}"
                    )
                else:
                    hover.append(
                        f"Protein: {item}<br>"
                        f"Sample: {group}<br>"
                        f"{state.metric_name}: {abs_value:{count_fmt}}<br>"
                        f"{_log_hover_line(state, abs_value)}"
                        f"Relative Contribution: {rel_value:.1f}%"
                    )

            all_values.extend(values)
            bar_color = color_map.get(bar_group, '#999')
            if bar_group == 'Other Proteins':
                bar_color = '#808080'
            # Replicate-based SEM error bars (absolute metrics only).
            bar_error_y = None
            if not state.is_relative and any(s > 0 for s in sems):
                bar_error_y = dict(type='data', array=sems, visible=True,
                                   thickness=1.2, width=3, color='#000000')
            fig.add_trace(go.Bar(
                x=x_pos, y=values, name=disp_bg,
                marker=dict(color=bar_color, line=dict(color='black', width=0.5)),
                width=bar_width * 0.9,
                error_y=bar_error_y,
                hovertext=hover, hoverinfo='text',
            ))

        if skipped_log_zero:
            state.warnings.append(
                f'{skipped_log_zero} bar(s) had a zero or missing value and were left blank '
                'under log transform (log of zero/negative is undefined). Turn off log '
                'transform to display them.'
            )
        y_title = state.ylabel or (f'Log<sub>{_log_sub(log_base)}</sub> ' if use_log else '') + state.metric_name
        finite_values = [v for v in all_values if v is not None]
        yaxis_prot = _build_grouped_yaxis(state, use_log, use_count, max(finite_values, default=0))
        x_lbl, legend_lbl = _orientation_axis_titles(state)
        fig.update_layout(
            **_common_layout(state),
            xaxis=dict(tickvals=list(range(len(categories))), ticktext=display_cats,
                       tickangle=45, title=x_lbl, **_axis_style(state.font_size_xlabel, state.font_size_xtick)),
            yaxis=dict(title=y_title, **yaxis_prot),
            barmode='group',
            hoverlabel=HOVERLABEL,
            legend=_legend_style(state, legend_lbl, yanchor='top', y=1.0, xanchor='left', x=1.05),
        )
        if show_sig and orientation == 'By Protein' and sig_geom:
            _add_grouped_significance(
                fig, categories, list(bar_groups_list), sig_geom,
                lambda name: state.protein_reps_dict.get(name, {}), metric_key,
                method=state.significance_method, sig_font_size=state.font_size_value_label)
        elif show_sig and orientation == 'By Sample' and bysample_tops:
            # Samples are the x-clusters; compare them on their replicate-level
            # filtered totals (total_reps_dict, computed on the selected-protein
            # filtered_df) and annotate per cluster. Numeric x-axis -> x_categorical=False.
            _add_totals_significance(fig, categories, state.total_reps_dict, metric_key,
                                     bysample_tops, is_log_axis=False,
                                     method=state.significance_method, x_categorical=False,
                                     sig_font_size=state.font_size_value_label)
        elif state.show_significance and state.is_relative:
            _add_significance_caption(fig, RELATIVE_STATS_NOTE)
        return fig

    elif plot_filter in ('No Filter', 'Both') and orientation == 'By Sample':
        # Sample totals: absolute bars, or per-sample % of the grand total when the
        # Relative metric is selected (plot_total_peptides handles both).
        return plot_total_peptides(state)

    return None


# ---------------------------------------------------------------------------
# Plot 3: Stacked bar plot
# ---------------------------------------------------------------------------

def _plot_no_filter_stacked(state: DataAnalysisState):
    """Stacked bar for 'No Filter': single 'Total' bar with one segment per sample group."""
    import plotly.graph_objects as go

    data = state.total_peptide_results_dict
    if not data:
        return None

    groups = list(data.keys())
    use_count = state.use_count
    color_seq = get_color_sequence(len(groups), state.color_scheme)

    total_count = sum(data[g]['unique_peptides'] for g in groups) or 1
    total_abund = sum(data[g]['total_Abundance'] for g in groups) or 1

    fig = go.Figure()
    for i, g in enumerate(groups):
        abs_val = data[g]['unique_peptides'] if use_count else data[g]['total_Abundance']
        rel_val = (data[g]['unique_peptides'] / total_count * 100) if use_count \
                  else (data[g]['total_Abundance'] / total_abund * 100)
        y_val = rel_val if state.is_relative else abs_val
        metric_label = 'Peptide Count' if use_count else 'Abundance'
        hover = (
            f"Sample: {g}<br>"
            f"{metric_label}: {abs_val:{state.num_format}}<br>"
            f"Contribution: {rel_val:.2f}%"
        )
        fig.add_trace(go.Bar(
            name=g, x=['Total'], y=[y_val],
            marker=dict(color=color_seq[i], line=dict(color='white', width=0.5)),
            hovertext=[hover], hoverinfo='text',
        ))

    y_title = ('Relative ' if state.is_relative else '') + state.metric_name
    yaxis_kw = _axis_style(state.font_size_ylabel, state.font_size_ytick)
    # Notebook parity (data_analysis.ipynb Plotter.plot_stacked_bar_scaled):
    # this bar stacks one segment per sample group, so absolute/scaled tick
    # values don't correspond to any single meaningful quantity — only the
    # relative (%) case gets readable tick labels.
    yaxis_kw['showticklabels'] = state.is_relative
    if state.is_relative:
        yaxis_kw['range'] = [0, 100]
        yaxis_kw['tickformat'] = '.1f'
    elif use_count:
        yaxis_kw['tickformat'] = ',d'
        yaxis_kw['exponentformat'] = 'none'
    else:
        # Abundance totals are large; use scientific notation ("3E+11") rather than
        # Plotly's default SI suffix ("300B"), matching the other stacked/bar plots.
        yaxis_kw['exponentformat'] = 'E'
        yaxis_kw['showexponent'] = 'all'
    fig.update_layout(
        **_common_layout(state),
        barmode='stack',
        xaxis=dict(title='', **_axis_style(state.font_size_xlabel, state.font_size_xtick)),
        yaxis=dict(title=state.ylabel or y_title, **yaxis_kw),
        hoverlabel=HOVERLABEL,
        legend=_legend_style(state, state.legend_title or 'Sample'),
    )
    return fig


def plot_stacked_bar_scaled(state: DataAnalysisState):
    """
    Stacked bar plot matching notebook's plot_stacked_bar_scaled logic.

    By Sample orientation (x=groups, legend=items):
      - y_value = rel_value / 100 * total_sums[g]  where total_sums uses the
        proper per-group total (unique functional peptide count / actual functional
        abundance) to prevent double-counting multi-function peptides.

    By Function / By Protein orientation (x=items, legend=groups):
      - Uses function_distribution_dict / protein_sample_distribution_dict with
        'relative' key = group's contribution to the item's total.
      - y_value = (rel_percentage / 100) * item_total  for non-log absolute
      - y_value = (rel_percentage / 100) * log10(item_total)  for log
    """
    import plotly.graph_objects as go

    selected_groups = state.selected_groups
    use_count = state.use_count
    use_log = state.log_transform
    log_base = state.log_base
    plot_filter = state.plot_filter
    orientation = state.orientation

    # "No Filter" is orientation-driven for stacked bars. By Sample keeps the single
    # stacked "Total" bar (sample composition). By Protein / By Function instead
    # plot EVERY protein / function — i.e. the identical figure to
    # "Selected Protein(s) + By Protein" / "Selected Function(s) + By Function" with
    # all items selected. (Previously every No-Filter case returned the Total bar,
    # so switching orientation appeared to do nothing.)
    all_items_override = None
    if plot_filter == 'No Filter':
        if orientation == 'By Protein':
            plot_filter = 'Selected Protein(s)'
            if state.protein_df is not None and not state.protein_df.empty:
                all_items_override = list(state.protein_df['Description'])
        elif orientation == 'By Function':
            plot_filter = 'Selected Function(s)'
            if state.function_df is not None and not state.function_df.empty:
                all_items_override = list(state.function_df['Description'])
        else:  # By Sample
            return _plot_no_filter_stacked(state)

    # ── Resolve data source and selected items ─────────────────────────────────
    is_func_filter = plot_filter in ('Selected Function(s)', 'Functional vs Non-Functional Peptides')
    is_prot_filter = plot_filter == 'Selected Protein(s)'
    is_both = plot_filter == 'Both'

    if is_func_filter:
        df = state.function_df
        items_col = 'Description'
        selected_items = all_items_override if all_items_override is not None else state.selected_functions
        item_label = 'Function'
    elif is_prot_filter:
        df = state.protein_df
        items_col = 'Description'
        selected_items = all_items_override if all_items_override is not None else state.selected_protein_names
        item_label = 'Protein'
    else:  # Both
        # Combine function and protein items
        fn_df = state.function_df
        pr_df = state.protein_df
        if (fn_df is None or fn_df.empty) and (pr_df is None or pr_df.empty):
            return plot_total_peptides(state)
        df = pd.concat([d for d in [fn_df, pr_df] if d is not None and not d.empty], ignore_index=True)
        items_col = 'Description'
        selected_items = list(state.selected_functions) + list(state.selected_protein_names)
        item_label = 'Item'

    if df is None or (hasattr(df, 'empty') and df.empty):
        return plot_total_peptides(state)

    # ── Build color map ────────────────────────────────────────────────────────
    color_seq = get_color_sequence(len(selected_items), state.color_scheme)
    color_map = {item: color_seq[i] for i, item in enumerate(selected_items)}
    if state.plot_minor and 'Other Functions' in color_map:
        color_map['Other Functions'] = '#808080'
    if state.plot_minor and 'Other Proteins' in color_map:
        color_map['Other Proteins'] = '#808080'

    fig = go.Figure()

    # ── By Sample orientation ──────────────────────────────────────────────────
    if orientation == 'By Sample':
        # Compute total_sums per group – matches notebook's approach of using the
        # real total (unique peptides / unique functional abundance) to prevent
        # stacked bars from summing to an inflated (double-counted) total.
        total_sums = {}
        for g in selected_groups:
            if use_count:
                if is_func_filter:
                    total_sums[g] = state.function_count_totals_dict.get(g, 0)
                elif is_prot_filter:
                    _pct = state.protein_count_bysample_dict.get(g, 0)
                    if not _pct and state.protein_df is not None and f'Count_{g}' in state.protein_df.columns:
                        _pct = int(state.protein_df[f'Count_{g}'].sum())
                    total_sums[g] = _pct
                else:  # Both
                    total_sums[g] = state.abundance_count_by_sample_dict.get(g, {}).get('unique_peptides', 0)
            else:  # abundance
                if is_func_filter:
                    total_sums[g] = state.function_abundance_totals_dict.get(g, 0)
                elif is_prot_filter:
                    sm = state.sum_df
                    col_key = f'Avg_{g}'
                    if sm is not None and col_key in sm['Sample'].values:
                        total_sums[g] = float(sm.loc[sm['Sample'] == col_key, 'Total_Sum'].values[0])
                    elif state.protein_df is not None and f'Avg_{g}' in state.protein_df.columns:
                        total_sums[g] = float(state.protein_df[f'Avg_{g}'].sum())
                    else:
                        total_sums[g] = 0.0
                else:  # Both
                    total_sums[g] = state.abundance_count_by_sample_dict.get(g, {}).get('total_Abundance', 0)

        for i, item in enumerate(selected_items):
            row = df[df[items_col] == item] if df is not None else pd.DataFrame()
            values = []
            hover_texts = []
            for g in selected_groups:
                if df is None or row.empty:
                    values.append(0)
                    hover_texts.append(f"No data for {item} in {g}")
                    continue

                abs_col = f'Count_{g}' if use_count else f'Avg_{g}'
                rel_col = f'Rel_Count_{g}' if use_count else f'Rel_Avg_{g}'
                abs_value = float(row[abs_col].values[0]) if abs_col in row.columns else 0.0
                rel_value = float(row[rel_col].values[0]) if rel_col in row.columns else 0.0

                total = total_sums.get(g, 0)
                if state.is_relative:
                    v = rel_value
                elif use_log:
                    log_total = _safe_log(total, log_base) if total > 0 else 0.0
                    v = rel_value / 100.0 * log_total
                else:
                    # Scale proportionally so bars sum to the true total
                    # (prevents double-counting of multi-function peptides)
                    v = rel_value / 100.0 * total if total > 0 else abs_value

                values.append(v)
                hover_texts.append(
                    f"{item_label}: {item}<br>"
                    f"Sample: {g}<br>"
                    f"Relative {state.metric_name}: {rel_value:.2f}%<br>"
                    f"Absolute {state.metric_name}: {abs_value:{state.num_format}}"
                )

            disp_name = redact_string_descriptions(item)
            color = color_map.get(item, '#999')
            if item in ('Other Functions', 'Other Proteins'):
                color = '#808080'
            fig.add_trace(go.Bar(
                name=disp_name, x=selected_groups, y=values,
                marker=dict(color=color, line=dict(color='white', width=0.5)),
                hovertext=hover_texts, hoverinfo='text',
            ))

        # Totals annotation (matches notebook's scatter-text trace above bars)
        if not state.is_relative:
            text_vals = []
            for g in selected_groups:
                total = total_sums.get(g, 0)
                if use_log and total > 0:
                    text_vals.append(f"{_log_arr(max(total, 1e-10), log_base):.2f}")
                elif use_count:
                    text_vals.append(f"{int(total):,}")
                else:
                    text_vals.append(f"{total:.2e}")
            y_totals = []
            for g in selected_groups:
                total = total_sums.get(g, 0)
                if use_log and total > 0:
                    y_totals.append(float(_log_arr(max(total, 1e-10), log_base)))
                else:
                    y_totals.append(total)
            fig.add_trace(go.Scatter(
                x=selected_groups, y=y_totals,
                mode='text', text=text_vals,
                textposition='top center',
                textfont=dict(size=max(8, state.font_size_value_label - 2) if len(selected_groups) > 12
                              else state.font_size_value_label, color='black'),
                showlegend=True,
                name=f'Total {state.metric_name}',
                hoverinfo='none',
            ))

    # ── By Function / By Protein orientation ──────────────────────────────────
    else:
        # Build per-item totals dict from function_distribution_dict or protein_sample_distribution_dict
        # Uses 'relative'  = each group's contribution to the item's total
        item_data_dict = {}
        for item in selected_items:
            if plot_filter in ('Selected Function(s)', 'Functional vs Non-Functional Peptides', 'Both'):
                d = state.function_distribution_dict.get(item)
                if d is not None:
                    item_data_dict[item] = d
            if item not in item_data_dict and plot_filter in ('Selected Protein(s)', 'Both'):
                d = state.protein_sample_distribution_dict.get(item)
                if d is None:
                    # try matching by name from protein_df
                    if state.protein_df is not None:
                        pr = state.protein_df[state.protein_df['Description'] == item]
                        if not pr.empty:
                            row = pr.iloc[0]
                            tot_ab = sum(float(row.get(f'Avg_{g}', 0)) for g in selected_groups)
                            tot_ct = sum(float(row.get(f'Count_{g}', 0)) for g in selected_groups)
                            d = {
                                'Abundance': {g: float(row.get(f'Avg_{g}', 0)) for g in selected_groups},
                                'counts': {g: float(row.get(f'Count_{g}', 0)) for g in selected_groups},
                                'relative': {
                                    g: float(row.get(f'Avg_{g}', 0)) / tot_ab * 100 if tot_ab > 0 else 0.0
                                    for g in selected_groups
                                },
                                'count_relative_to_function': {
                                    g: float(row.get(f'Count_{g}', 0)) / tot_ct * 100 if tot_ct > 0 else 0.0
                                    for g in selected_groups
                                },
                                'total_Abundance': tot_ab,
                                'unique_peptide_count': int(row.get('unique_peptide_count', 0)),
                            }
                if d is not None:
                    item_data_dict[item] = d

        all_items = [it for it in selected_items if it in item_data_dict]
        display_items = [redact_string_descriptions(it) for it in all_items]
        colors = get_color_sequence(len(selected_groups), state.color_scheme)

        for gi, g in enumerate(selected_groups):
            values = []
            hover_texts = []
            for item in all_items:
                data = item_data_dict[item]

                # Use group's contribution to item total ('relative' key)
                # so stacked bars show distribution of this item across samples.
                # For counts the per-group share lives under a source-specific
                # key: 'count_relative_to_function' for functions,
                # 'count_relative_to_protein' for proteins (and a "Both" figure
                # mixes the two). Fall back across both so protein bars don't
                # collapse to zero. See protein_sample_distribution_dict /
                # function_distribution_dict in data_processor.
                if use_count:
                    rel = data.get('count_relative_to_function')
                    if rel is None:
                        rel = data.get('count_relative_to_protein', {})
                else:
                    rel = data.get('relative', {})
                rel_percentage = rel.get(g, 0.0)

                abs_value = (
                    data.get('counts' if use_count else 'Abundance', {}).get(g, 0.0)
                )
                item_total = (
                    data.get('unique_peptide_count', 0)
                    if use_count
                    else data.get('total_Abundance', 0.0)
                )

                if state.is_relative:
                    v = rel_percentage
                elif use_log:
                    log_total = _safe_log(item_total, log_base) if item_total > 0 else 0.0
                    v = rel_percentage / 100.0 * log_total
                else:
                    # Equivalent to abs_value but expressed via proportional total
                    # (same as notebook: total_Abundance * rel_percentage / 100)
                    v = item_total * (rel_percentage / 100.0) if item_total > 0 else abs_value

                values.append(v)
                abs_count_label = 'Count' if use_count else 'Abundance'
                hover_texts.append(
                    f"{item_label}: {item}<br>"
                    f"Sample: {g}<br>"
                    f"Sample's contribution: {rel_percentage:.2f}%<br>"
                    f"{abs_count_label} in sample: {abs_value:{state.num_format}}"
                )

            color = colors[gi] if gi < len(colors) else '#CCCCCC'
            fig.add_trace(go.Bar(
                name=redact_string_descriptions(g),
                x=display_items,
                y=values,
                marker=dict(color=color, line=dict(color='white', width=0.5)),
                hovertext=hover_texts, hoverinfo='text',
            ))

        # Totals above bars
        if not state.is_relative:
            text_vals = []
            y_totals = []
            for item in all_items:
                data = item_data_dict[item]
                total = (
                    data.get('unique_peptide_count', 0)
                    if use_count
                    else data.get('total_Abundance', 0.0)
                )
                if use_log and total > 0:
                    text_vals.append(f"{_log_arr(max(total, 1e-10), log_base):.2f}")
                    y_totals.append(float(_log_arr(max(total, 1e-10), log_base)))
                elif use_count:
                    text_vals.append(f"{int(total):,}")
                    y_totals.append(total)
                else:
                    text_vals.append(f"{total:.2e}")
                    y_totals.append(total)
            fig.add_trace(go.Scatter(
                x=display_items, y=y_totals,
                mode='text', text=text_vals,
                textposition='top center',
                textfont=dict(size=state.font_size_value_label, color='black'),
                showlegend=True,
                name=f'Total {state.metric_name}',
                hoverinfo='none',
            ))

    # ── Layout ─────────────────────────────────────────────────────────────────
    y_title = ('Relative ' if state.is_relative else '') + \
               (f'Log<sub>{_log_sub(log_base)}</sub> ' if (use_log and not state.is_relative) else '') + state.metric_name
    yaxis_stk = _axis_style(state.font_size_ylabel, state.font_size_ytick)
    # Notebook parity (data_analysis.ipynb Plotter.plot_stacked_bar_scaled): once
    # items are stacked into layered segments, absolute/scaled tick values don't
    # map to any one meaningful number — only the relative (%) case shows ticks.
    yaxis_stk['showticklabels'] = state.is_relative
    if state.is_relative:
        yaxis_stk['range'] = [0, 100]
        yaxis_stk['tickformat'] = '.1f'
    elif use_count and not use_log:
        yaxis_stk['tickformat'] = ',d'
        yaxis_stk['exponentformat'] = 'none'
        yaxis_stk['showexponent'] = 'none'
    elif use_count and use_log:
        yaxis_stk['tickformat'] = '.1f'
        yaxis_stk['exponentformat'] = 'none'
    elif use_log:
        yaxis_stk['tickformat'] = '.1f'
        yaxis_stk['exponentformat'] = 'none'
    else:
        # Abundance, absolute, non-log: the stacked total (already log-transformed
        # abundance data included) can still be small, in which case plain integers
        # read better than scientific notation. See LINEAR_SCALE_THRESHOLD.
        if orientation == 'By Sample':
            stack_max = max(total_sums.values(), default=0)
        else:
            stack_max = max((d.get('total_Abundance', 0.0) for d in item_data_dict.values()), default=0)
        if stack_max < LINEAR_SCALE_THRESHOLD:
            yaxis_stk['tickformat'] = ',d'
            yaxis_stk['exponentformat'] = 'none'
            yaxis_stk['showexponent'] = 'none'
        else:
            yaxis_stk['exponentformat'] = 'E'
            yaxis_stk['showexponent'] = 'all'

    xlabel, legend_lbl = _orientation_axis_titles(state)
    fig.update_layout(**{
        **_common_layout(state),
        'barmode': 'stack',
        'xaxis': dict(title=xlabel, tickangle=-90 if orientation == 'By Sample' else 45,
                      **_axis_style(state.font_size_xlabel, state.font_size_xtick)),
        'yaxis': dict(title=state.ylabel or y_title, **yaxis_stk),
        'hoverlabel': HOVERLABEL,
        'legend': _legend_style(state, legend_lbl, yanchor='top', y=0.95, xanchor='left', x=1.05,
                                bgcolor='rgba(255,255,255,0.9)'),
        'height': 820,
        'margin': dict(t=100, l=100, r=100, b=100),
    })
    return fig


# ---------------------------------------------------------------------------
# Plot 4: Pie charts
# ---------------------------------------------------------------------------

def create_pie_charts(state: DataAnalysisState):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    selected_groups = state.selected_groups
    use_count = state.use_count
    plot_filter = state.plot_filter
    orientation = state.orientation

    if plot_filter in ('Selected Function(s)', 'Functional vs Non-Functional Peptides'):
        df = state.function_df
        items_col = 'Description'
        items = state.selected_functions
    elif plot_filter == 'Selected Protein(s)':
        df = state.protein_df
        items_col = 'Description'
        items = state.selected_protein_names
    else:
        # No Filter / Both → one pie per sample showing total distribution
        metrics = state.abundance_count_by_sample_dict
        labels = list(metrics.keys())
        vals = [v['unique_peptides'] if use_count else v['total_Abundance'] for v in metrics.values()]
        colors = get_color_sequence(len(labels), state.color_scheme)
        fig = go.Figure()
        _hover_fmt_nf = "Sample: %{label}<br>Value: %{value:,.0f}<br><extra></extra>" if use_count else "Sample: %{label}<br>Value: %{value:.2e}<br><extra></extra>"
        fig.add_trace(go.Pie(
            labels=labels, values=vals,
            marker_colors=colors,
            textposition='inside', textinfo='percent',
            textfont=dict(size=state.font_size_value_label),
            hovertemplate=_hover_fmt_nf,
        ))
        fig.update_layout(**{**_common_layout(state),
                          'showlegend': True,
                          'hoverlabel': HOVERLABEL,
                          'legend': _legend_style(state, state.legend_title or 'Sample')})
        return fig

    if df is None or df.empty:
        return None

    if orientation == 'By Sample':
        # One pie per sample
        n = len(selected_groups)
        cols = min(n, 3)
        rows = (n + cols - 1) // cols
        fig = make_subplots(rows=rows, cols=cols, specs=[[{'type': 'pie'}] * cols] * rows)

        for i, group in enumerate(selected_groups):
            row_idx = i // cols + 1
            col_idx = i % cols + 1
            col_key = (f'Count_{group}' if use_count else f'Avg_{group}')
            if col_key in df.columns:
                pie_df = df[df[items_col].isin(items) & (df[col_key] > 0)]
            else:
                pie_df = df[df[items_col].isin(items)]
            labels = list(pie_df[items_col])
            vals = [float(v) for v in pie_df.get(col_key, pd.Series([0] * len(pie_df)))]
            colors = get_color_sequence(len(labels), state.color_scheme)
            colors = [
                '#808080' if lbl in ('Other Proteins', 'Other Functions') else c
                for lbl, c in zip(labels, colors)
            ]

            hover_fmt = "%{label}: %{value:,.0f}<br><extra></extra>" if use_count else "%{label}: %{value:.2e}<br><extra></extra>"
            fig.add_trace(go.Pie(
                labels=[redact_string_descriptions(l) for l in labels],
                values=vals, name=group,
                marker_colors=colors,
                textposition='inside', textinfo='percent',
                textfont=dict(size=state.font_size_value_label),
                title=dict(text=group, font=dict(size=14)),
                hovertemplate=hover_fmt,
            ), row=row_idx, col=col_idx)

        fig.update_layout(**{**_common_layout(state),
                          'height': 400 * rows, 'width': 1000,
                          'hoverlabel': HOVERLABEL,
                          'legend': _legend_style(state, state.legend_title or 'Item')})
        return fig
    else:
        # By Function/Protein: one pie per item
        items_to_plot = [it for it in items if it in df[items_col].values]
        n = len(items_to_plot)
        if n == 0:
            return None
        cols = min(n, 3)
        rows = (n + cols - 1) // cols
        fig = make_subplots(rows=rows, cols=cols, specs=[[{'type': 'pie'}] * cols] * rows)

        for i, item in enumerate(items_to_plot):
            row_idx = i // cols + 1
            col_idx = i % cols + 1
            labels = selected_groups
            row = df[df[items_col] == item]
            vals = [float(row.get(f'Count_{g}', 0).values[0] if not row.empty else 0)
                    if use_count else
                    float(row.get(f'Avg_{g}', pd.Series([0])).values[0] if not row.empty else 0)
                    for g in selected_groups]
            colors = get_color_sequence(len(labels), state.color_scheme)

            _item_hover_fmt = "%{label}: %{value:,.0f}<br><extra></extra>" if use_count else "%{label}: %{value:.2e}<br><extra></extra>"
            fig.add_trace(go.Pie(
                labels=labels, values=vals, name=item,
                marker_colors=colors,
                textposition='inside', textinfo='percent',
                textfont=dict(size=state.font_size_value_label),
                title=dict(text=redact_string_descriptions(item), font=dict(size=12)),
                hovertemplate=_item_hover_fmt,
            ), row=row_idx, col=col_idx)

        fig.update_layout(**{**_common_layout(state),
                          'height': 400 * rows, 'width': 1000,
                          'hoverlabel': HOVERLABEL,
                          'legend': _legend_style(state, state.legend_title or 'Sample')})
        return fig


# ---------------------------------------------------------------------------
# Plot 5: Correlation scatter
# ---------------------------------------------------------------------------

def create_correlation_plot(state: DataAnalysisState):
    from scipy.stats import pearsonr, spearmanr
    import plotly.graph_objects as go

    selected_groups = state.selected_groups
    df = state.filtered_df
    if df is None or len(selected_groups) != 2:
        return None

    g1, g2 = selected_groups[0], selected_groups[1]
    col1, col2 = f'Avg_{g1}', f'Avg_{g2}'
    if col1 not in df.columns or col2 not in df.columns:
        return None

    use_log = state.log_transform
    first_color = get_single_color(state.color_scheme)

    fdf = df[(df[col1] > 0) & (df[col2] > 0)].copy()
    if fdf.empty:
        return None

    if use_log:
        base = state.log_base
        x_vals = _log_arr(fdf[col1], base)
        y_vals = _log_arr(fdf[col2], base)
        tickfmt = '.0f'
        sub = _log_sub(base)
        x_label, y_label = f'Log<sub>{sub}</sub> ({g1})', f'Log<sub>{sub}</sub> ({g2})'
    else:
        x_vals = fdf[col1]
        y_vals = fdf[col2]
        # Already-log-transformed abundance data can still land here with small
        # values; scientific notation is misleading noise below the same
        # linear-scale threshold used elsewhere. See LINEAR_SCALE_THRESHOLD.
        if max(float(x_vals.max()), float(y_vals.max())) < LINEAR_SCALE_THRESHOLD:
            tickfmt = ',d'
        else:
            tickfmt = '.1e'
        x_label, y_label = g1, g2

    corr_text = 'n/a'
    if len(fdf) > 1:
        if state.correlation_type == 'Pearson':
            corr, p = pearsonr(x_vals, y_vals)
            sym = '<i>r</i>'
        else:
            corr, p = spearmanr(x_vals, y_vals)
            sym = 'ρ'
        p_str = 'p < 0.001' if p < 0.001 else f'p = {p:.3g}'
        corr_text = f'{sym} = {corr:.3f}, {p_str}'

    # Hover columns
    id_col = 'Unique Peptide ID'
    fn_col = 'function' if 'function' in fdf.columns else None
    pn_col = 'protein_name' if 'protein_name' in fdf.columns else None

    customdata = np.column_stack([
        fdf[id_col] if id_col in fdf.columns else [''] * len(fdf),
        fdf[fn_col].fillna('N/A') if fn_col else ['N/A'] * len(fdf),
        fdf[pn_col].fillna('N/A') if pn_col else ['N/A'] * len(fdf),
    ])

    # Synchronized tick values for both axes (same range)
    x_min, x_max = float(x_vals.min()), float(x_vals.max())
    y_min, y_max = float(y_vals.min()), float(y_vals.max())
    overall_min = min(x_min, y_min)
    overall_max = max(x_max, y_max)
    combined_range = [overall_min * 0.95, overall_max * 1.05]
    if use_log:
        range_span = combined_range[1] - combined_range[0]
        if range_span > 5:
            tick_vals = list(np.arange(np.floor(combined_range[0]), np.ceil(combined_range[1]) + 1, 1))
        else:
            tick_vals = list(np.arange(np.floor(combined_range[0] * 2) / 2,
                                       np.ceil(combined_range[1] * 2) / 2 + 0.5, 0.5))
    else:
        tick_vals = list(np.linspace(combined_range[0], combined_range[1], 6))

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=x_vals, y=y_vals, mode='markers',
        name=f'Correlation: {corr_text}',
        marker=dict(color=first_color),
        customdata=customdata,
        hovertemplate=(
            '<b>Peptide ID:</b> %{customdata[0]}<br>'
            '<b>Function:</b> %{customdata[1]}<br>'
            '<b>Protein:</b> %{customdata[2]}<br>'
            '<b>%{xaxis.title.text}:</b> %{x:' + tickfmt + '}<br>'
            '<b>%{yaxis.title.text}:</b> %{y:' + tickfmt + '}<br>'
            '<extra></extra>'
        ),
    ))

    if len(fdf) > 1:
        z = np.polyfit(x_vals, y_vals, 1)
        x_range_vals = np.linspace(x_vals.min(), x_vals.max(), 100)
        fig.add_trace(go.Scatter(
            x=x_range_vals, y=np.poly1d(z)(x_range_vals),
            mode='lines', line=dict(color=first_color, dash='dash'),
            name='Trendline', showlegend=True,
            hovertemplate='<extra></extra>',
        ))

    ax_kw = dict(tickformat=tickfmt, tickvals=tick_vals, nticks=5,
                 range=combined_range,
                 showline=True, linewidth=1, linecolor='black',
                 gridcolor='lightgray', showgrid=True, zeroline=False)

    leg_title = state.legend_title or f'{state.correlation_type} Correlation Values'
    fig.update_layout(
        title=dict(text=_make_title(state, kind='correlation'), x=0.5, xanchor='center',
                   font=dict(size=state.font_size_plot_title, color='black')),
        xaxis=dict(title=x_label, tickangle=45,
                   title_font=dict(size=state.font_size_xlabel, color='black'),
                   tickfont=dict(size=state.font_size_xtick, color='black'), **ax_kw),
        yaxis=dict(title=y_label,
                   title_font=dict(size=state.font_size_ylabel, color='black'),
                   tickfont=dict(size=state.font_size_ytick, color='black'), **ax_kw),
        height=500, width=600,
        template='plotly_white',
        font=BASE_FONT,
        hoverlabel=HOVERLABEL,
        showlegend=True,
        legend=_legend_style(state, leg_title, yanchor='top', y=0.99, xanchor='right', x=0.99,
                             bgcolor='rgba(255, 255, 255, 0.8)'),
        margin=dict(t=100, b=80, l=80, r=50),
    )
    return fig


def create_correlation_splom(state: DataAnalysisState):
    from scipy.stats import pearsonr, spearmanr
    import plotly.graph_objects as go

    selected_groups = state.selected_groups
    df = state.filtered_df
    if df is None or len(selected_groups) < 3:
        return None

    use_log = state.log_transform
    base = state.log_base
    sub = _log_sub(base)
    first_color = get_single_color(state.color_scheme)
    ct = state.correlation_type

    # Filter to positive values in all groups
    fdf = df.copy()
    valid_groups = []
    for g in selected_groups:
        col = f'Avg_{g}'
        if col in fdf.columns:
            fdf = fdf[fdf[col] > 0]
            valid_groups.append(g)

    if fdf.empty or len(valid_groups) < 3:
        return None

    dimensions = []
    all_vals = []
    for g in valid_groups:
        col = f'Avg_{g}'
        vals = _log_arr(fdf[col], base) if use_log else fdf[col]
        dimensions.append(dict(values=vals, label=f'Log<sub>{sub}</sub> ({g})' if use_log else g))
        all_vals.extend(vals)

    id_col = 'Unique Peptide ID'
    fn_col = 'function' if 'function' in fdf.columns else None
    pn_col = 'protein_name' if 'protein_name' in fdf.columns else None
    customdata = np.column_stack([
        fdf[id_col] if id_col in fdf.columns else [''] * len(fdf),
        fdf[fn_col].fillna('N/A') if fn_col else ['N/A'] * len(fdf),
        fdf[pn_col].fillna('N/A') if pn_col else ['N/A'] * len(fdf),
    ])

    if use_log:
        tickfmt = '.0f'
    else:
        # Already-log-transformed abundance data can still land here with small
        # values; scientific notation is misleading noise below the same
        # linear-scale threshold used elsewhere. See LINEAR_SCALE_THRESHOLD.
        tickfmt = ',d' if max(all_vals) < LINEAR_SCALE_THRESHOLD else '.1e'
    splom = go.Splom(
        dimensions=dimensions,
        marker=dict(color=first_color, size=8, line=dict(width=1, color='white')),
        diagonal=dict(visible=False),
        showupperhalf=False,
        customdata=customdata,
        hovertemplate=(
            '<b>Peptide:</b> %{customdata[0]}<br>'
            '<b>Function:</b> %{customdata[1]}<br>'
            '<b>Protein:</b> %{customdata[2]}<br>'
            '<extra></extra>'
        ),
        showlegend=False,
    )

    corr_traces = []
    for i, g1 in enumerate(valid_groups):
        for j, g2 in enumerate(valid_groups):
            if i >= j:
                continue
            c1, c2 = f'Avg_{g1}', f'Avg_{g2}'
            x_vals = _log_arr(fdf[c1], base) if use_log else fdf[c1]
            y_vals = _log_arr(fdf[c2], base) if use_log else fdf[c2]
            if ct == 'Pearson':
                corr, p = pearsonr(x_vals, y_vals)
                sym = 'r'
            else:
                corr, p = spearmanr(x_vals, y_vals)
                sym = 'ρ'
            p_str = 'p < 0.001' if p < 0.001 else f'p = {p:.3g}'
            corr_traces.append(go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(color=first_color),
                name=f'{g1}–{g2}: {sym} = {corr:.3f}, {p_str}',
                showlegend=True,
            ))

    fig = go.Figure(data=[splom] + corr_traces)
    # Cap the overall figure so many groups shrink the panels rather than
    # growing the figure past the viewport. Uncapped, this scaled linearly
    # (20 groups -> 5000px) and overflowed the browser window.
    sz = min(250 * len(valid_groups), MAX_SPLOM_SIZE)
    leg_title = state.legend_title or f'{ct} Correlation'
    fig.update_layout(
        title=dict(text=_make_title(state, kind='correlation'), x=0.5, xanchor='center',
                   font=dict(size=state.font_size_plot_title, color='black')),
        width=sz + 250, height=sz,
        template='plotly_white',
        font=BASE_FONT,
        hoverlabel=HOVERLABEL,
        legend=_legend_style(state, leg_title, x=1.01, y=1, yanchor='top'),
    )
    # Apply consistent axis ranges/formats
    overall_min, overall_max = min(all_vals), max(all_vals)
    combined_range = [overall_min * 0.95, overall_max * 1.05]
    tick_vals = np.linspace(combined_range[0], combined_range[1], 6)
    for k in range(1, len(valid_groups) + 1):
        fig.update_layout({
            f'xaxis{k}': dict(range=combined_range, tickvals=tick_vals, tickformat=tickfmt,
                              tickfont=dict(size=state.font_size_xtick, color='black'), tickangle=45,
                              title_font=dict(size=state.font_size_xlabel, color='black')),
            f'yaxis{k}': dict(range=combined_range, tickvals=tick_vals, tickformat=tickfmt,
                              tickfont=dict(size=state.font_size_ytick, color='black'),
                              title_font=dict(size=state.font_size_ylabel, color='black')),
        })
    return fig


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_grouped_yaxis(state: DataAnalysisState, use_log: bool, use_count: bool, max_value: float = None) -> dict:
    """Return y-axis kwargs matching notebook grouped bar plot formatting."""
    kw = _axis_style(state.font_size_ylabel, state.font_size_ytick)
    if state.is_relative:
        kw['range'] = [0, 100]
    elif use_count and use_log:
        kw['type'] = 'linear'
        kw['tickformat'] = '.1f'
        kw['exponentformat'] = 'none'
        kw['showexponent'] = 'none'
    elif use_count and not use_log:
        kw['type'] = 'linear'
        kw['tickformat'] = ',d'
        kw['exponentformat'] = 'none'
        kw['showexponent'] = 'none'
    elif not use_count and use_log:
        kw['type'] = 'linear'
        kw['tickformat'] = '.1f'
        kw['exponentformat'] = 'none'
        kw['showexponent'] = 'none'
    elif max_value is not None and max_value < LINEAR_SCALE_THRESHOLD:
        # Small values (including abundance data already log-transformed
        # upstream) read better as plain integers than scientific notation.
        kw['type'] = 'linear'
        kw['tickformat'] = ',d'
        kw['exponentformat'] = 'none'
        kw['showexponent'] = 'none'
    else:  # abundance, no log
        kw['type'] = 'linear'
        kw['exponentformat'] = 'E'
        kw['showexponent'] = 'all'
    return kw


def _power_ticks(fig, y_axis='yaxis', x_axis=None):
    """
    Replace numeric tick labels on log-transformed axes with 10^n HTML superscript format.
    Scans all trace y (and optionally x) values to determine integer tick range.
    """
    import math

    def _vals_from_traces(attr):
        out = []
        for trace in fig.data:
            raw = getattr(trace, attr, None)
            if raw is None:
                continue
            try:
                out.extend(v for v in raw if v is not None and not (isinstance(v, float) and math.isnan(v)))
            except TypeError:
                pass
        return out

    def _make_ticks(vals):
        if not vals:
            return None, None
        mn, mx = math.floor(min(vals)), math.ceil(max(vals))
        tv = list(range(mn, mx + 1))
        tt = [f'10<sup>{v}</sup>' if v >= 0 else f'10<sup>−{abs(v)}</sup>' for v in tv]
        return tv, tt

    y_vals = _vals_from_traces('y')
    tv, tt = _make_ticks(y_vals)
    if tv:
        fig.update_layout({y_axis: dict(tickvals=tv, ticktext=tt)})

    if x_axis:
        x_vals = _vals_from_traces('x')
        tv, tt = _make_ticks(x_vals)
        if tv:
            fig.update_layout({x_axis: dict(tickvals=tv, ticktext=tt)})


def _common_layout(state: DataAnalysisState) -> dict:
    # NB: callers spread this as **kwargs and add their own hoverlabel/legend,
    # so those keys must NOT be set here (would raise a duplicate-kwarg error).
    return dict(
        template='plotly_white',
        # Width intentionally omitted: with config.responsive=True, Plotly.js
        # only stretches a dimension left unset in the layout. A fixed width
        # here pins the plot to that pixel size regardless of container size,
        # which is what produced the large blank strip on wide containers
        # (the plot rendered at 1000px inside a much wider #plot-container).
        # Height stays fixed since chart height isn't container-width-dependent.
        height=800,
        margin=PLOT_MARGIN,
        font=BASE_FONT,
        title=dict(text=_make_title(state), y=0.95, x=0.5,
                   xanchor='center', yanchor='top',
                   font=dict(size=state.font_size_plot_title, color='black')),
    )


def _axis_style(title_size: int = 18, tick_size: int = 16) -> dict:
    return dict(
        showline=True, linewidth=1, linecolor='black',
        mirror=False, gridcolor='lightgray', showgrid=True, zeroline=False,
        title_font=dict(size=title_size, color='black'),
        tickfont=dict(size=tick_size, color='black'),
    )


def _legend_style(state: DataAnalysisState, text: str, **extra) -> dict:
    """Legend dict with the Appearance Settings legend-title font size driving
    both the legend title and the legend item text, per the Appearance Settings
    control (one input governs the whole legend's type size)."""
    size = state.font_size_legend_title
    return dict(
        title=dict(text=text, font=dict(size=size, color='black')),
        font=dict(size=size, color='black'),
        **extra,
    )


def _orientation_axis_titles(state: DataAnalysisState):
    """Return (x_axis_label, legend_title) for the current orientation, matching
    the notebook's invert_plot / x-axis-label logic. When the x-axis holds the
    categories (By Function / By Protein) it is named for them and the legend
    holds the samples; By Sample keeps samples on the x-axis and names the legend
    for whatever is stacked/grouped. User-supplied labels always win.
    """
    orient = state.orientation
    pf = state.plot_filter
    if orient == 'By Protein':
        x_lbl, legend = 'Proteins', 'Samples'
    elif orient == 'By Function':
        x_lbl, legend = 'Functions', 'Samples'
    else:  # By Sample
        x_lbl = 'Samples'
        if pf in ('Selected Function(s)', 'Functional vs Non-Functional Peptides'):
            legend = 'Functions'
        elif pf == 'Selected Protein(s)':
            legend = 'Proteins'
        else:
            legend = 'Sample'
    return (state.xlabel or x_lbl, state.legend_title or legend)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def generate_plot(merged_df: pd.DataFrame, group_data_dict: dict, protein_dict: dict, params: dict):
    """
    Main entry point: build state, run pipeline, generate appropriate plot.
    Returns (plotly_json_str, warnings_list).
    """
    warnings = []

    state = DataAnalysisState(merged_df, group_data_dict, protein_dict, params)

    try:
        state.run_pipeline()
    except Exception:
        warnings.append('Error during data processing: ' + traceback.format_exc())
        return None, warnings

    plot_type = state.plot_type
    orientation = state.orientation
    plot_filter = state.plot_filter
    use_count = state.use_count
    fig = None

    try:
        if plot_type == 'Grouped Bar Plots':
            if orientation == 'By Sample' and plot_filter in ('No Filter', 'Both') and not state.is_relative:
                fig = plot_total_peptides(state)
            else:
                fig = create_grouped_bar_plot(state)

        elif plot_type == 'Stacked Bar Plots':
            fig = plot_stacked_bar_scaled(state)

        elif plot_type == 'Pie Charts':
            fig = create_pie_charts(state)

        elif plot_type == 'Corr. Scatter Plots':
            n = len(state.selected_groups)
            if n < 2:
                warnings.append('Please select at least 2 groups for correlation analysis.')
            elif n == 2:
                fig = create_correlation_plot(state)
            else:
                fig = create_correlation_splom(state)
    except Exception:
        warnings.append('Error generating plot: ' + traceback.format_exc())
        return None, warnings

    warnings.extend(state.warnings)

    if fig is None:
        warnings.append('No plot generated. Please check your selections and data.')
        return None, warnings

    return fig.to_json(), warnings
