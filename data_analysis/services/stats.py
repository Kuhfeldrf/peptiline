"""Group-comparison statistics for the Data Analysis bar charts.

Compares the selected sample groups on the *replicate-level* summed value of
whatever quantity the bar shows (summed abundance or unique-peptide count), so
the statistic tests exactly the quantity plotted rather than pseudo-replicating
across peptides. For two groups the result is a significance bracket
(``*``/``**``/``***``); for three or more it is a compact-letter display (groups
that share a letter are not significantly different).

Two methods are offered (``method=`` on :func:`compare_groups`):

* ``'tukey'`` (default) — one-way **ANOVA + Tukey HSD**. This is the convention
  used across the Dallas-lab / Kuhfeld peptidomics papers, so the app reports
  comparisons the same way the authors' own manuscripts do. It assumes roughly
  equal within-group variance.
* ``'games-howell'`` — **Welch's ANOVA + Games–Howell** post-hoc, which does
  *not* assume equal variance. Peptidomic abundance data are typically
  heteroscedastic (higher-abundance groups vary more), and simulation at n≈3
  shows plain Tukey's false-positive rate inflating to ~2× nominal under a 5:1
  SD ratio while Games–Howell holds at ~0.05 — at a modest, conservative power
  cost when variances happen to be equal. Prefer this when group spreads differ.

Kraskal–Wallis / Dunn was evaluated and rejected for this setting: at n≈3 its
false-positive rate is poorly controlled (~0.10–0.22 in simulation), it tests
stochastic dominance rather than a location shift, and it yields coarse discrete
p-values — a worse fit than the two variance-based options above.

**Power is the binding constraint, not the choice of test.** At the replicate
counts typical here (n≈3) even a large true effect (Cohen's d≈1) is detected
only ~15% of the time; you need d≈3 for ~80% power. A shared letter / ``ns``
therefore means "not resolvable at this sample size", NOT "no difference".
A comparison is only reported when every compared group has at least
``MIN_REPLICATES`` replicates; under-replicated groups are *excluded and named*
(see the ``excluded`` key) rather than silently dropped.
"""
import warnings
from itertools import combinations

import numpy as np
from scipy.stats import f, f_oneway, studentized_range, tukey_hsd

# Minimum replicates per group required before any significance is reported.
# Below this a t-test/ANOVA is too underpowered to be meaningful at these sample
# sizes, so we deliberately show effect size (bars + SEM) without p-values.
MIN_REPLICATES = 3
ALPHA = 0.05

# Accepted method aliases -> canonical key.
_METHOD_ALIASES = {
    'tukey': 'tukey', 'anova': 'tukey', 'anova_tukey': 'tukey',
    'anova+tukey': 'tukey', 'anova + tukey hsd': 'tukey',
    'games-howell': 'games-howell', 'games_howell': 'games-howell',
    'gameshowell': 'games-howell', 'welch': 'games-howell',
    'welch-gh': 'games-howell', 'welch + games-howell': 'games-howell',
}
_METHOD_LABELS = {
    'tukey': 'ANOVA + Tukey HSD',
    'games-howell': 'Welch ANOVA + Games–Howell',
}


def significance_stars(p, alpha=ALPHA):
    """Map a p-value to a significance label (``*`` / ``**`` / ``***`` / ``ns``)."""
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return 'ns'
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < alpha:
        return '*'
    return 'ns'


def _bron_kerbosch(adj, nodes):
    """All maximal cliques of an undirected graph (adjacency as {node: set}).

    Nodes here are group indices and an edge means "not significantly
    different". Each maximal clique becomes one letter in the compact-letter
    display: two groups share a letter iff they are adjacent (not different),
    which is exactly the property a correct CLD needs. Graphs are tiny (a handful
    of groups), so the plain recursive algorithm is more than adequate.
    """
    cliques = []

    def expand(r, p, x):
        if not p and not x:
            cliques.append(r)
            return
        for v in list(p):
            expand(r | {v}, p & adj[v], x & adj[v])
            p = p - {v}
            x = x | {v}

    expand(set(), set(nodes), set())
    return cliques


def _compact_letters(groups, pairwise, alpha=ALPHA):
    """Compact-letter display for a set of groups given pairwise Tukey p-values.

    ``groups`` is an ordered list; ``pairwise`` maps frozenset({g1, g2}) -> p.
    Returns {group: 'a'} where groups sharing any letter are NOT significantly
    different.
    """
    idx = {g: i for i, g in enumerate(groups)}
    k = len(groups)
    # adjacency = "not significantly different" (self-loops excluded)
    adj = {i: set() for i in range(k)}
    for i in range(k):
        for j in range(i + 1, k):
            p = pairwise.get(frozenset((groups[i], groups[j])))
            not_diff = not (p is not None and not np.isnan(p) and p < alpha)
            if not_diff:
                adj[i].add(j)
                adj[j].add(i)

    cliques = _bron_kerbosch(adj, range(k))
    # Deterministic letter order: sort cliques by their smallest member.
    cliques.sort(key=lambda c: min(c))
    letters = {g: '' for g in groups}
    for c, clique in enumerate(cliques):
        ch = chr(ord('a') + c) if c < 26 else f'({c + 1})'
        for i in clique:
            letters[groups[i]] += ch
    return letters


def letters_are_uninformative(letters):
    """True when every group carries the *same* compact letter.

    A CLD in which all groups share one identical letter says only "nothing was
    resolved" — but drawn on the figure it reads as a per-bar annotation and
    invites the reader to hunt for a contrast that does not exist (the exact
    confusion reported for the multi-function bar chart). Callers use this to
    drop the letter row and state the result in the caption instead.
    """
    values = [v for v in letters.values() if v]
    return bool(values) and len(values) == len(letters) and len(set(values)) == 1


def _welch_anova(arrays):
    """Welch's (unequal-variance) one-way ANOVA omnibus p-value.

    Returns NaN if any group has zero within-group variance (the weights n/var
    are then undefined). The omnibus is informational only — the markers are
    driven by the pairwise Games–Howell p-values, which stay well-defined for a
    constant group as long as its partner varies.
    """
    k = len(arrays)
    n = np.array([len(a) for a in arrays], dtype=float)
    m = np.array([a.mean() for a in arrays])
    v = np.array([a.var(ddof=1) for a in arrays])
    if np.any(v == 0) or np.any(~np.isfinite(v)):
        return float('nan')
    w = n / v
    W = w.sum()
    mbar = (w * m).sum() / W
    num = (w * (m - mbar) ** 2).sum() / (k - 1)
    tmp = ((1 - w / W) ** 2 / (n - 1)).sum()
    denom = 1 + (2 * (k - 2) / (k ** 2 - 1)) * tmp
    F = num / denom
    df2 = (k ** 2 - 1) / (3 * tmp)
    return float(f.sf(F, k - 1, df2))


def _games_howell(arrays, groups):
    """Pairwise Games–Howell p-values: {frozenset({g_i, g_j}): p}.

    The Welch-style post-hoc for Tukey — uses per-pair Welch–Satterthwaite
    degrees of freedom against the studentised-range distribution, so it does not
    assume equal variance or equal n. A pair where one group is constant is still
    well defined (its variance contributes 0) provided the other group varies; a
    pair where *both* are constant yields NaN (handled upstream by the
    zero-within-variance guard).
    """
    k = len(arrays)
    m = [a.mean() for a in arrays]
    v = [a.var(ddof=1) for a in arrays]
    n = [len(a) for a in arrays]
    out = {}
    for i, j in combinations(range(k), 2):
        vni, vnj = v[i] / n[i], v[j] / n[j]
        se = np.sqrt(vni + vnj)
        key = frozenset((groups[i], groups[j]))
        if se == 0 or not np.isfinite(se):
            out[key] = float('nan')
            continue
        q = abs(m[i] - m[j]) / se * np.sqrt(2)
        df = (vni + vnj) ** 2 / (vni ** 2 / (n[i] - 1) + vnj ** 2 / (n[j] - 1))
        out[key] = float(studentized_range.sf(q, k, df))
    return out


def compare_groups(rep_arrays, alpha=ALPHA, method='tukey'):
    """Compare sample groups on their replicate-level values.

    Parameters
    ----------
    rep_arrays : dict
        ``{group_name: [replicate values]}`` in display order.
    alpha : float
        Significance threshold for the compact-letter display.
    method : str
        ``'tukey'`` (default; ANOVA + Tukey HSD, assumes equal variance) or
        ``'games-howell'`` (Welch ANOVA + Games–Howell, robust to unequal
        variance). Aliases accepted (see ``_METHOD_ALIASES``); anything
        unrecognised falls back to ``'tukey'``.

    Returns
    -------
    dict with keys:
        testable : bool      – whether a valid comparison was produced
        reason   : str|None  – why not, when ``testable`` is False
        groups   : list      – groups that entered the test (n >= MIN_REPLICATES)
        excluded : list      – groups dropped for too few replicates (named, not
                               silently discarded)
        n        : dict       – replicate count per group (after dropping NaNs)
        min_n    : int|None   – smallest replicate count among tested groups
                               (small values => low power; a caption cue)
        anova_p  : float|None – omnibus p-value
        pairwise : dict       – {frozenset({g1, g2}): p}
        letters  : dict       – {group: 'a'} compact-letter display (>= 3 groups)
        method   : str        – human-readable description for the figure caption
    """
    method = _METHOD_ALIASES.get(str(method).strip().lower(), 'tukey')
    ordered = list(rep_arrays.keys())
    # Drop NaN/inf replicate values BEFORE counting n, so a group with a missing
    # replicate is counted at its real size and never yields a NaN p-value that
    # masquerades as a valid ("testable") result.
    cleaned = {}
    for g in ordered:
        a = np.asarray(rep_arrays[g], dtype=float)
        cleaned[g] = a[np.isfinite(a)]
    n = {g: int(cleaned[g].size) for g in ordered}
    base = {
        'testable': False, 'reason': None, 'groups': [], 'excluded': [],
        'n': n, 'min_n': None, 'anova_p': None, 'pairwise': {}, 'letters': {},
        'method': _METHOD_LABELS[method],
    }

    valid = [g for g in ordered if n[g] >= MIN_REPLICATES]
    # Groups present but under-replicated: reported so the UI can say what was
    # left out rather than the comparison quietly ignoring them.
    base['excluded'] = [g for g in ordered if 0 < n[g] < MIN_REPLICATES]
    if len(valid) < 2:
        base['reason'] = f'need ≥{MIN_REPLICATES} replicates in ≥2 groups'
        return base

    arrays = [cleaned[g] for g in valid]
    # Degenerate: every value identical across all groups -> no test possible.
    if np.ptp(np.concatenate(arrays)) == 0:
        base['reason'] = 'no variance among replicates'
        return base

    # When every group is internally constant (zero within-group variance) the
    # pooled standard error is 0, so both ANOVA/Tukey and Welch/Games–Howell are
    # undefined (÷0). scipy does NOT return NaN for Tukey here — f_oneway gives
    # p=0.0 and tukey_hsd divides by zero to an off-diagonal p of 0.0 — which
    # would fabricate a "significant" result from a degenerate input. Guard it
    # explicitly so we never claim a significance we cannot legitimately compute.
    if max(float(np.var(a, ddof=1)) for a in arrays) == 0:
        base['reason'] = 'zero within-group variance (test undefined)'
        return base

    # Any remaining ÷0 / invalid warnings from partially-degenerate inputs are
    # suppressed; the resulting NaN p-values fall through to 'ns'.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with np.errstate(divide='ignore', invalid='ignore'):
            if method == 'games-howell':
                try:
                    anova_p = _welch_anova(arrays)
                except Exception:
                    anova_p = None
                try:
                    pairwise = _games_howell(arrays, valid)
                except Exception:
                    pairwise = {}
            else:
                try:
                    anova_p = float(f_oneway(*arrays).pvalue)
                except Exception:
                    anova_p = None
                pairwise = {}
                try:
                    pv = tukey_hsd(*arrays).pvalue
                    for i, j in combinations(range(len(valid)), 2):
                        pairwise[frozenset((valid[i], valid[j]))] = float(pv[i, j])
                except Exception:
                    pairwise = {}

    letters = _compact_letters(valid, pairwise, alpha) if len(valid) >= 3 else {}
    return {
        'testable': True, 'reason': None, 'groups': valid,
        'excluded': base['excluded'], 'n': n,
        'min_n': min(n[g] for g in valid),
        'anova_p': anova_p, 'pairwise': pairwise, 'letters': letters,
        'method': _METHOD_LABELS[method],
    }


# ---------------------------------------------------------------------------
# Two-group effect sizes for the heatmap differential-comparison track
# (SMD / log2FC). These are point estimates over already-aggregated replicate
# vectors — the significance machinery above is for the bar charts; here we
# report effect magnitude per position/peptide. See
# manuscript/docs/HEATMAP_DIFFERENTIAL_STATS.md.
# ---------------------------------------------------------------------------

# Minimum finite replicates per group for a defined Cohen's d. A pooled SD needs
# a within-group variance, so each side needs >= 2 values. This is the maths
# floor for the descriptive track, deliberately below MIN_REPLICATES (which
# gates *significance*, not effect size).
MIN_REP_FOR_D = 2

# Sawilowsky (2009) "new effect size rules of thumb", descending. Used only for
# the human-readable magnitude label in captions/legends.
_SAWILOWSKY = (
    (2.0, 'huge'), (1.2, 'very large'), (0.8, 'large'),
    (0.5, 'medium'), (0.2, 'small'), (0.01, 'very small'),
)


def cohens_d(a, b):
    """Cohen's *d* standardized mean difference between two replicate samples.

    ``d = (mean_a - mean_b) / s_pooled`` with the classic (n-1)-weighted pooled
    standard deviation. **Positive d => higher in ``a``** (the first-selected
    series). Non-finite values are dropped listwise before counting n / mean /
    SD (the "drop" missing-value policy).

    Returns ``nan`` — so callers draw a gap rather than a spurious value — when
    either group has fewer than ``MIN_REP_FOR_D`` finite values or the pooled SD
    is 0 / non-finite (effect size undefined).
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    na, nb = a.size, b.size
    if na < MIN_REP_FOR_D or nb < MIN_REP_FOR_D:
        return float('nan')
    va, vb = a.var(ddof=1), b.var(ddof=1)
    s_pooled = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    if s_pooled == 0 or not np.isfinite(s_pooled):
        return float('nan')
    return float((a.mean() - b.mean()) / s_pooled)


def log2_fold_change(a, b, eps=0.0):
    """log2 fold change of group means: ``log2((mean_a + eps)/(mean_b + eps))``.

    Positive => higher in ``a``. Non-finite values dropped listwise. ``eps`` is a
    pseudocount to keep the ratio defined when a mean is 0; with ``eps=0`` any
    zero/negative mean yields ``nan`` (rendered as a gap). This is a magnitude of
    change, NOT a variance-aware effect size — pair it with :func:`cohens_d`.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return float('nan')
    ma, mb = a.mean() + eps, b.mean() + eps
    if ma <= 0 or mb <= 0 or not np.isfinite(ma) or not np.isfinite(mb):
        return float('nan')
    return float(np.log2(ma / mb))


def effect_size_label(d):
    """Map \\|Cohen's d\\| to the Sawilowsky (2009) magnitude label."""
    if d is None or (isinstance(d, float) and np.isnan(d)):
        return 'n/a'
    ad = abs(float(d))
    for threshold, label in _SAWILOWSKY:
        if ad >= threshold:
            return label
    return 'negligible'
