"""
Unit tests for the Data Analysis group-comparison statistics
(peptide.data_analysis.services.stats): ANOVA + Tukey HSD, significance stars,
and the compact-letter display used to annotate the bar charts.

The final class is an end-to-end check that the significance overlay actually
renders onto the bar charts (processor -> plotter) when ``show_significance`` is
set, so the wiring between the stats service and the figures stays covered.
"""
import unittest

import numpy as np
import pandas as pd

# conftest.py bootstraps Django; import the service directly.
from data_analysis.services import stats
from data_analysis.services import plotter
from data_analysis.services.data_processor import DataAnalysisState


class TestSignificanceStars(unittest.TestCase):
    def test_thresholds(self):
        self.assertEqual(stats.significance_stars(0.0005), '***')
        self.assertEqual(stats.significance_stars(0.005), '**')
        self.assertEqual(stats.significance_stars(0.03), '*')
        self.assertEqual(stats.significance_stars(0.2), 'ns')

    def test_none_and_nan(self):
        self.assertEqual(stats.significance_stars(None), 'ns')
        self.assertEqual(stats.significance_stars(float('nan')), 'ns')


class TestCompareGroups(unittest.TestCase):
    def test_two_groups_clearly_different(self):
        res = stats.compare_groups({'A': [10, 11, 9, 10], 'B': [100, 102, 98, 101]})
        self.assertTrue(res['testable'])
        self.assertLess(res['anova_p'], 0.05)
        p = res['pairwise'][frozenset(('A', 'B'))]
        self.assertEqual(stats.significance_stars(p), '***')

    def test_two_groups_not_different(self):
        res = stats.compare_groups({'A': [10, 11, 9, 10], 'B': [10, 9, 11, 10]})
        self.assertTrue(res['testable'])
        self.assertEqual(stats.significance_stars(res['pairwise'][frozenset(('A', 'B'))]), 'ns')

    def test_insufficient_replicates(self):
        res = stats.compare_groups({'A': [10, 11], 'B': [100, 102]})
        self.assertFalse(res['testable'])
        self.assertIn('replicates', res['reason'])

    def test_no_variance(self):
        res = stats.compare_groups({'A': [5, 5, 5], 'B': [5, 5, 5]})
        self.assertFalse(res['testable'])

    def test_zero_within_group_variance_is_safe(self):
        # Each group internally constant but different -> pooled standard error is
        # 0, so ANOVA/Tukey are undefined (÷0). scipy returns p=0.0 here rather
        # than NaN, so we must guard explicitly: not testable, and no fabricated
        # significance. Must not raise and must not over-claim.
        res = stats.compare_groups({'A': [5, 5, 5], 'B': [9, 9, 9]})
        self.assertFalse(res['testable'])
        self.assertIn('variance', res['reason'])
        for p in res['pairwise'].values():
            self.assertEqual(stats.significance_stars(p), 'ns')

    def test_one_constant_group_still_tests(self):
        # Only ONE group internally constant -> pooled variance > 0, so the test
        # is well defined and must still run (the zero-variance guard is scoped to
        # the all-groups-constant case, not a single flat group).
        res = stats.compare_groups({'A': [5, 5, 5], 'B': [9, 10, 11]})
        self.assertTrue(res['testable'])
        self.assertEqual(
            stats.significance_stars(res['pairwise'][frozenset(('A', 'B'))]), '***')


class TestNaNAndExclusion(unittest.TestCase):
    """Missing values and under-replicated groups must be handled, not hidden."""

    def test_nan_replicate_dropped_and_recounted(self):
        # A has a NaN -> only 2 real values -> below MIN_REPLICATES. It must be
        # dropped and counted at its true size, NOT enter the test at full n and
        # yield a NaN p-value masquerading as a valid comparison.
        res = stats.compare_groups({'A': [10, float('nan'), 9], 'B': [100, 99, 101]})
        self.assertEqual(res['n']['A'], 2)
        self.assertIn('A', res['excluded'])
        self.assertFalse(res['testable'])  # only B left with >=3 -> <2 valid groups
        # No fabricated significance survives.
        for p in res['pairwise'].values():
            self.assertEqual(stats.significance_stars(p), 'ns')

    def test_underreplicated_group_is_named_not_silently_dropped(self):
        res = stats.compare_groups({'A': [10, 11], 'B': [50, 51, 49], 'C': [90, 91, 89]})
        self.assertTrue(res['testable'])
        self.assertEqual(res['groups'], ['B', 'C'])
        self.assertEqual(res['excluded'], ['A'])   # surfaced for the caption

    def test_min_n_tracks_smallest_tested_group(self):
        res = stats.compare_groups({'A': [1, 2, 3], 'B': [10, 11, 12, 13, 14]})
        self.assertEqual(res['min_n'], 3)


class TestGamesHowell(unittest.TestCase):
    """Welch ANOVA + Games–Howell: the unequal-variance alternative."""

    def test_method_label_and_aliases(self):
        self.assertEqual(
            stats.compare_groups({'A': [1, 2, 3], 'B': [4, 5, 6]}, method='games-howell')['method'],
            'Welch ANOVA + Games–Howell')
        # aliases resolve to the same method
        for alias in ('welch', 'games_howell', 'Welch + Games-Howell'):
            self.assertEqual(
                stats.compare_groups({'A': [1, 2, 3]}, method=alias)['method'],
                'Welch ANOVA + Games–Howell')
        # unknown -> default tukey
        self.assertEqual(
            stats.compare_groups({'A': [1, 2, 3]}, method='bogus')['method'],
            'ANOVA + Tukey HSD')

    def test_clear_difference_detected(self):
        res = stats.compare_groups({'A': [10, 11, 9], 'B': [100, 99, 101]},
                                   method='games-howell')
        self.assertTrue(res['testable'])
        self.assertEqual(
            stats.significance_stars(res['pairwise'][frozenset(('A', 'B'))]), '***')

    def test_no_difference_not_flagged(self):
        res = stats.compare_groups({'A': [10, 11, 9, 10], 'B': [10, 9, 11, 10]},
                                   method='games-howell')
        self.assertEqual(
            stats.significance_stars(res['pairwise'][frozenset(('A', 'B'))]), 'ns')

    def test_compact_letters_three_groups(self):
        res = stats.compare_groups(
            {'A': [10, 11, 9, 10], 'B': [11, 10, 12, 9], 'C': [200, 205, 195, 201]},
            method='games-howell')
        letters = res['letters']
        self.assertTrue(set(letters['A']) & set(letters['B']))   # A,B share
        self.assertFalse(set(letters['A']) & set(letters['C']))  # A,C differ

    def test_controls_type_i_error_under_heteroscedasticity(self):
        # The whole reason this method exists: under strongly unequal variance and
        # equal means, plain Tukey over-calls significance (~2x nominal at n=3)
        # while Games–Howell stays near alpha. Verify GH's empirical false-positive
        # rate is well controlled and materially lower than Tukey's here.
        rng = np.random.default_rng(2024)
        n_sim, n = 400, 3
        fp_tukey = fp_gh = 0
        for _ in range(n_sim):
            data = {'A': rng.normal(100, 2, n), 'B': rng.normal(100, 30, n)}
            pt = stats.compare_groups(data, method='tukey')['pairwise'].get(frozenset(('A', 'B')))
            pg = stats.compare_groups(data, method='games-howell')['pairwise'].get(frozenset(('A', 'B')))
            fp_tukey += stats.significance_stars(pt) != 'ns'
            fp_gh += stats.significance_stars(pg) != 'ns'
        rate_tukey, rate_gh = fp_tukey / n_sim, fp_gh / n_sim
        self.assertLess(rate_gh, 0.10, f'GH FPR too high: {rate_gh:.3f}')
        self.assertLess(rate_gh, rate_tukey + 1e-9,
                        f'GH ({rate_gh:.3f}) should not exceed Tukey ({rate_tukey:.3f})')


class TestCompactLetters(unittest.TestCase):
    def test_two_clusters(self):
        # A≈B, C far apart -> A,B share a letter; C stands alone.
        res = stats.compare_groups({
            'A': [10, 11, 9, 10], 'B': [11, 10, 12, 9], 'C': [200, 205, 195, 201],
        })
        letters = res['letters']
        self.assertTrue(set(letters['A']) & set(letters['B']))         # A, B share
        self.assertFalse(set(letters['A']) & set(letters['C']))        # A, C differ
        self.assertFalse(set(letters['B']) & set(letters['C']))        # B, C differ

    def test_all_different(self):
        res = stats.compare_groups({
            'A': [10, 11, 9], 'B': [60, 62, 58], 'C': [200, 205, 195],
        })
        letters = res['letters']
        # No two groups share any letter.
        for g1, g2 in (('A', 'B'), ('A', 'C'), ('B', 'C')):
            self.assertFalse(set(letters[g1]) & set(letters[g2]), f'{g1},{g2} should not share')

    def test_all_same(self):
        res = stats.compare_groups({
            'A': [10, 11, 9], 'B': [10, 12, 9], 'C': [11, 10, 10],
        })
        letters = res['letters']
        # All mutually non-different -> all share one common letter.
        common = set(letters['A']) & set(letters['B']) & set(letters['C'])
        self.assertTrue(common)


class TestSignificanceOverlayIntegration(unittest.TestCase):
    """Processor -> plotter: the significance overlay renders on the figures."""

    GROUPS = ['Ctrl', 'LowDose', 'HighDose']

    @staticmethod
    def _letters(fig):
        return [a for a in fig.layout.annotations if a.text and '<b>' in a.text]

    @staticmethod
    def _caption(fig):
        # The significance footnote is appended to the x-axis title (below the
        # tick labels via automargin), tagged with a grey 11px span.
        try:
            title = fig.layout.xaxis.title.text or ''
        except Exception:
            title = ''
        return title if 'color:#555555' in title else None

    def _dataset(self, nreps=3):
        rng = np.random.default_rng(0)
        reps = {g: [f'{g}_{i}' for i in range(1, nreps + 1)] for g in self.GROUPS}
        rows = []
        for i in range(12):
            fn = 'ACE-inhibitory' if i % 2 == 0 else 'Antioxidant'
            prot = 'P02666' if i < 6 else 'P02662'
            row = {'Unique Peptide ID': f'pep{i}', 'Protein': prot, 'function': fn}
            for g in self.GROUPS:
                base = {'Ctrl': 100, 'LowDose': 120, 'HighDose': 500}[g]
                for rc in reps[g]:
                    row[rc] = float(base + rng.normal(0, 8))
                row[f'Avg_{g}'] = float(np.mean([row[rc] for rc in reps[g]]))
            rows.append(row)
        protein_dict = {'P02666': {'name': 'Beta-casein'},
                        'P02662': {'name': 'Alpha-S1-casein'}}
        return pd.DataFrame(rows), reps, protein_dict

    def _state(self, reps, merged, protein_dict, **extra):
        params = dict(selected_groups=self.GROUPS,
                      selected_functions=['All Functional Peptides'],
                      selected_proteins=['All Proteins (No Filter)'],
                      abs_or_count='Abundance', metric_type='Absolute',
                      show_significance=True)
        params.update(extra)
        st = DataAnalysisState(merged, reps, protein_dict, params)
        st.run_pipeline()
        return st

    def _has_error_bars(self, fig):
        return any(t.error_y is not None and t.error_y.array is not None for t in fig.data)

    def test_grouped_bars_get_letters_error_bars_and_caption(self):
        merged, reps, pdct = self._dataset()
        for orient, pfilter, log in [
            ('By Function', 'Selected Function(s)', False),
            ('By Function', 'Selected Function(s)', True),   # log axis path
            ('By Protein', 'Selected Protein(s)', False),
        ]:
            st = self._state(reps, merged, pdct, orientation=orient,
                             plot_filter=pfilter, log_transform=log)
            fig = plotter.create_grouped_bar_plot(st)
            # 2 categories x 3 groups -> one compact letter annotation per bar.
            self.assertEqual(len(self._letters(fig)), 6,
                             f'{orient} log={log}: expected 6 CLD letters')
            self.assertTrue(self._has_error_bars(fig),
                            f'{orient} log={log}: expected replicate SEM error bars')
            cap = self._caption(fig)
            self.assertIsNotNone(cap, f'{orient}: expected a method/power caption')
            self.assertIn('Tukey', cap)                  # default method
            self.assertIn('n = 3', cap)                  # low-power note at n=3

    def test_games_howell_method_selected_end_to_end(self):
        merged, reps, pdct = self._dataset()
        st = self._state(reps, merged, pdct, orientation='By Function',
                         plot_filter='Selected Function(s)',
                         significance_method='games-howell')
        fig = plotter.create_grouped_bar_plot(st)
        self.assertEqual(len(self._letters(fig)), 6)
        self.assertIn('Games', self._caption(fig))

    def test_caption_explains_when_not_enough_replicates(self):
        # Only 2 replicates per group -> nothing testable. The overlay must not be
        # silent: a caption explains WHY no significance is shown.
        merged, reps, pdct = self._dataset(nreps=2)
        st = self._state(reps, merged, pdct, orientation='By Function',
                         plot_filter='Selected Function(s)')
        fig = plotter.create_grouped_bar_plot(st)
        self.assertEqual(len(self._letters(fig)), 0)
        cap = self._caption(fig)
        self.assertIsNotNone(cap)
        self.assertIn('not shown', cap)
        self.assertIn('replicate', cap)

    def test_relative_metric_suppresses_overlay_and_error_bars(self):
        # A percentage bar's replicate SEM is not the SEM of the ratio, and the
        # significance overlay is defined on absolute values, so Relative mode must
        # carry neither error bars nor letters nor a caption.
        merged, reps, pdct = self._dataset()
        st = self._state(reps, merged, pdct, orientation='By Function',
                         plot_filter='Selected Function(s)', metric_type='Relative')
        fig = plotter.create_grouped_bar_plot(st)
        self.assertEqual(len(fig.layout.annotations), 0)
        self.assertFalse(self._has_error_bars(fig))

    def test_totals_plot_gets_compact_letters(self):
        merged, reps, pdct = self._dataset()
        st = self._state(reps, merged, pdct, orientation='By Sample',
                         plot_filter='No Filter')
        fig = plotter.plot_total_peptides(st)
        # One text trace for the numeric labels, one for the CLD letters.
        text_traces = [t for t in fig.data if getattr(t, 'mode', None) == 'text']
        self.assertGreaterEqual(len(text_traces), 2)
        self.assertIn('Tukey', self._caption(fig))


class TestAxisAndScaleFeatures(unittest.TestCase):
    """Tasks 2 & 4: linear-below-10^4 default and orientation-aware axis titles."""

    GROUPS = ['Ctrl', 'LowDose', 'HighDose']

    def _dataset(self, base_scale=1.0):
        rng = np.random.default_rng(1)
        reps = {g: [f'{g}_{i}' for i in (1, 2, 3)] for g in self.GROUPS}
        rows = []
        for i in range(12):
            fn = 'ACE-inhibitory' if i % 2 == 0 else 'Antioxidant'
            prot = 'P02666' if i < 6 else 'P02662'
            row = {'Unique Peptide ID': f'pep{i}', 'Protein': prot, 'function': fn}
            for g in self.GROUPS:
                base = {'Ctrl': 100, 'LowDose': 120, 'HighDose': 150}[g] * base_scale
                for rc in reps[g]:
                    row[rc] = float(base + rng.normal(0, 5))
                row[f'Avg_{g}'] = float(np.mean([row[rc] for rc in reps[g]]))
            rows.append(row)
        return (pd.DataFrame(rows), reps,
                {'P02666': {'name': 'Beta-casein'}, 'P02662': {'name': 'Alpha-S1-casein'}})

    def _state(self, base_scale=1.0, **extra):
        merged, reps, pdct = self._dataset(base_scale)
        params = dict(selected_groups=self.GROUPS,
                      selected_functions=['All Functional Peptides'],
                      selected_proteins=['All Proteins (No Filter)'],
                      abs_or_count='Abundance', metric_type='Absolute',
                      plot_filter='No Filter', orientation='By Sample',
                      log_transform=False)
        params.update(extra)
        st = DataAnalysisState(merged, reps, pdct, params)
        st.run_pipeline()
        return st

    def test_task2_linear_axis_below_threshold(self):
        # Per-group summed abundance here is ~1.2k-1.8k (< 10^4) and log is off,
        # so the totals abundance axis must default to linear, not log.
        st = self._state(base_scale=1.0)
        fig = plotter.plot_total_peptides(st)
        self.assertEqual(fig.layout.yaxis.type, 'linear')

    def test_task2_log_axis_above_threshold(self):
        # Scale abundance well above 10^4 -> keep the default log axis.
        st = self._state(base_scale=100.0)   # ~120k-180k per group
        fig = plotter.plot_total_peptides(st)
        self.assertEqual(fig.layout.yaxis.type, 'log')

    def test_task2_log_checkbox_still_works_below_threshold(self):
        # The checkbox is independent of the auto-default: ticking it log-transforms
        # even when values are small (task 2 must not disable the checkbox).
        st = self._state(base_scale=1.0, log_transform=True)
        fig = plotter.plot_total_peptides(st)
        # Values are pre-log-transformed onto a linear axis; the title says Log10.
        self.assertIn('Log', fig.layout.yaxis.title.text)

    def test_log2_transform_labels_and_scales_axis(self):
        # Log2 is mutually exclusive with Log10; when selected the axis title
        # carries the base-2 subscript and the plotted values are log2, not log10.
        st10 = self._state(base_scale=100.0, log_transform=True, log_base=10)
        st2 = self._state(base_scale=100.0, log_transform=True, log_base=2)
        fig10 = plotter.plot_total_peptides(st10)
        fig2 = plotter.plot_total_peptides(st2)
        self.assertIn('Log<sub>2</sub>', fig2.layout.yaxis.title.text)
        self.assertIn('Log<sub>10</sub>', fig10.layout.yaxis.title.text)
        y10 = max(v for v in fig10.data[0].y if v is not None)
        y2 = max(v for v in fig2.data[0].y if v is not None)
        # log2(x) == log10(x) / log10(2)  ->  ratio ~3.3219
        self.assertAlmostEqual(y2 / y10, 1.0 / np.log10(2), places=3)

    def test_task4_orientation_axis_titles(self):
        cases = {
            ('By Function', 'Selected Function(s)'): ('Functions', 'Samples'),
            ('By Protein', 'Selected Protein(s)'): ('Proteins', 'Samples'),
            ('By Sample', 'Selected Function(s)'): ('Samples', 'Functions'),
            ('By Sample', 'Selected Protein(s)'): ('Samples', 'Proteins'),
        }
        for (orient, pf), (want_x, want_leg) in cases.items():
            st = self._state(orientation=orient, plot_filter=pf)
            self.assertEqual(plotter._orientation_axis_titles(st), (want_x, want_leg),
                             f'{orient}/{pf}')

    def test_task4_user_labels_override(self):
        st = self._state(orientation='By Function', plot_filter='Selected Function(s)',
                         xlabel='Custom X', legend_title='Custom L')
        self.assertEqual(plotter._orientation_axis_titles(st), ('Custom X', 'Custom L'))


class TestSignificanceBracketPlacement(unittest.TestCase):
    """Task 5: brackets/letters sit above the bars, error bars, and value labels."""

    def _two_group_totals_fig(self):
        rng = np.random.default_rng(3)
        groups = ['A', 'B']
        reps = {g: [f'{g}_{i}' for i in (1, 2, 3)] for g in groups}
        rows = []
        for i in range(10):
            row = {'Unique Peptide ID': f'p{i}', 'Protein': 'P1', 'function': 'ACE-inhibitory'}
            for g in groups:
                for rc in reps[g]:
                    row[rc] = float(100 + rng.normal(0, 6))
                row[f'Avg_{g}'] = float(np.mean([row[rc] for rc in reps[g]]))
            rows.append(row)
        params = dict(selected_groups=groups, plot_filter='No Filter',
                      orientation='By Sample', abs_or_count='Abundance',
                      metric_type='Absolute', show_significance=True, log_transform=False)
        st = DataAnalysisState(pd.DataFrame(rows), reps, {'P1': {'name': 'P1'}}, params)
        st.run_pipeline()
        return plotter.plot_total_peptides(st), st

    def test_totals_bracket_is_above_bars_and_labels(self):
        fig, st = self._two_group_totals_fig()
        # Bracket = 3 line shapes; the crossbar y must clear the tallest bar+SEM
        # (the value labels sit right at bar+SEM, so above that clears them too).
        line_shapes = [s for s in fig.layout.shapes if s.type == 'line']
        self.assertGreaterEqual(len(line_shapes), 3)
        bar = next(t for t in fig.data if t.type == 'bar')
        sem = bar.error_y.array if bar.error_y and bar.error_y.array is not None else [0] * len(bar.y)
        max_bar_top = max(y + s for y, s in zip(bar.y, sem))
        crossbar_y = max(s.y0 for s in line_shapes)
        self.assertGreater(crossbar_y, max_bar_top,
                           'significance bracket must sit above the bar+error-bar top')
        # No stray numeric categories were introduced on the x-axis.
        self.assertEqual(list(bar.x), ['A', 'B'])

    def test_totals_no_numeric_categories_leaked(self):
        # Regression guard: earlier a scatter with numeric x added '0','1','0.5'
        # categories. Every x value across traces must be a real group name.
        fig, st = self._two_group_totals_fig()
        allowed = set(st.selected_groups)
        for t in fig.data:
            xs = getattr(t, 'x', None) or []
            for x in xs:
                self.assertIn(x, allowed, f'unexpected x category {x!r} from trace {t.type}')


class TestNoFilterOrientation(unittest.TestCase):
    """No-Filter stacked bars are orientation-driven (all proteins / all functions)."""

    GROUPS = ['Ctrl', 'LowDose', 'HighDose']

    def _state(self, **extra):
        rng = np.random.default_rng(5)
        reps = {g: [f'{g}_{i}' for i in (1, 2, 3)] for g in self.GROUPS}
        rows = []
        for i in range(12):
            fn = 'ACE-inhibitory' if i % 2 == 0 else 'Antioxidant'
            prot = ['P1', 'P2', 'P3'][i % 3]
            row = {'Unique Peptide ID': f'pep{i}', 'Protein': prot, 'function': fn}
            for g in self.GROUPS:
                for rc in reps[g]:
                    row[rc] = float(100 + rng.normal(0, 8))
                row[f'Avg_{g}'] = float(np.mean([row[rc] for rc in reps[g]]))
            rows.append(row)
        pdct = {'P1': {'name': 'Prot-1'}, 'P2': {'name': 'Prot-2'}, 'P3': {'name': 'Prot-3'}}
        params = dict(selected_groups=self.GROUPS, plot_type='Stacked Bar Plots',
                      selected_functions=['All Functional Peptides'],
                      selected_proteins=['All Proteins (No Filter)'],
                      abs_or_count='Abundance', metric_type='Absolute',
                      plot_filter='No Filter', log_transform=False)
        params.update(extra)
        st = DataAnalysisState(pd.DataFrame(rows), reps, pdct, params)
        st.run_pipeline()
        return st

    def test_by_protein_plots_all_proteins(self):
        st = self._state(orientation='By Protein')
        fig = plotter.plot_stacked_bar_scaled(st)
        self.assertEqual(fig.layout.barmode, 'stack')
        self.assertEqual(fig.layout.xaxis.title.text, 'Proteins')
        # Every protein present should be an x-category (not a single "Total" bar).
        n_prot = st.protein_df['Description'].nunique()
        self.assertGreaterEqual(n_prot, 3)
        bar = next(t for t in fig.data if t.type == 'bar')
        self.assertEqual(len(bar.x), n_prot)

    def test_by_function_plots_all_functions(self):
        st = self._state(orientation='By Function')
        fig = plotter.plot_stacked_bar_scaled(st)
        self.assertEqual(fig.layout.xaxis.title.text, 'Functions')
        bar = next(t for t in fig.data if t.type == 'bar')
        self.assertGreaterEqual(len(bar.x), 2)

    def test_by_sample_still_single_total_bar(self):
        st = self._state(orientation='By Sample')
        fig = plotter.plot_stacked_bar_scaled(st)
        bar = next(t for t in fig.data if t.type == 'bar')
        self.assertEqual(list(bar.x), ['Total'])

    def test_by_protein_count_bars_are_nonzero(self):
        # Regression: By Protein + Peptide Count zeroed every stacked segment
        # because the plotter looked up 'count_relative_to_function' (a function
        # key) on protein data keyed by 'count_relative_to_protein'. The bars
        # collapsed to 0 while the floating "Total" labels remained.
        st = self._state(orientation='By Protein', abs_or_count='Count')
        fig = plotter.plot_stacked_bar_scaled(st)
        bar_traces = [t for t in fig.data if t.type == 'bar']
        # At least one segment per protein must carry a nonzero count.
        self.assertTrue(any(any(v > 0 for v in t.y) for t in bar_traces))
        # Stacked segments must sum to each protein's total-count label.
        totals = next(t for t in fig.data
                      if t.type == 'scatter' and 'Total' in (t.name or ''))
        for col, expected in enumerate(totals.y):
            stacked = sum(t.y[col] for t in bar_traces)
            self.assertAlmostEqual(stacked, expected, places=6)


class TestRelativeMetricNoFilter(unittest.TestCase):
    """Regression: the Relative toggle must apply under No Filter / Both too."""

    GROUPS = ['Ctrl', 'LowDose', 'HighDose']

    def _state(self, **extra):
        rng = np.random.default_rng(7)
        reps = {g: [f'{g}_{i}' for i in (1, 2, 3)] for g in self.GROUPS}
        rows = []
        for i in range(12):
            row = {'Unique Peptide ID': f'pep{i}', 'Protein': 'P1',
                   'function': 'ACE-inhibitory'}
            for g in self.GROUPS:
                base = {'Ctrl': 1000, 'LowDose': 2000, 'HighDose': 3000}[g]
                for rc in reps[g]:
                    row[rc] = float(base + rng.normal(0, 50))
                row[f'Avg_{g}'] = float(np.mean([row[rc] for rc in reps[g]]))
            rows.append(row)
        params = dict(selected_groups=self.GROUPS, plot_type='Grouped Bar Plots',
                      selected_functions=['All Functional Peptides'],
                      selected_proteins=['All Proteins (No Filter)'],
                      abs_or_count='Abundance', metric_type='Relative',
                      plot_filter='No Filter', orientation='By Sample', log_transform=True)
        params.update(extra)
        st = DataAnalysisState(pd.DataFrame(rows), reps, {'P1': {'name': 'P1'}}, params)
        st.run_pipeline()
        return st

    def _bar(self, fig):
        return next(t for t in fig.data if t.type == 'bar')

    def test_no_filter_by_sample_relative_scales_to_percent(self):
        st = self._state(plot_filter='No Filter', metric_type='Relative')
        fig = plotter.plot_total_peptides(st)
        self.assertEqual(list(fig.layout.yaxis.range), [0, 100])
        ys = list(self._bar(fig).y)
        self.assertTrue(all(0 <= y <= 100 for y in ys), ys)
        self.assertAlmostEqual(sum(ys), 100.0, places=3)  # composition sums to 100%

    def test_both_by_sample_relative_scales_to_percent(self):
        st = self._state(plot_filter='Both', metric_type='Relative')
        fig = plotter.plot_total_peptides(st)
        self.assertEqual(list(fig.layout.yaxis.range), [0, 100])
        self.assertAlmostEqual(sum(self._bar(fig).y), 100.0, places=3)

    def test_no_filter_by_sample_absolute_unchanged(self):
        # Absolute still renders the totals chart (not percentage-scaled).
        st = self._state(plot_filter='No Filter', metric_type='Absolute', log_transform=False)
        fig = plotter.plot_total_peptides(st)
        ys = list(self._bar(fig).y)
        self.assertGreater(max(ys), 100)  # real abundance values, not percentages

    def test_relative_count_no_filter(self):
        st = self._state(plot_filter='No Filter', metric_type='Relative', abs_or_count='Count')
        fig = plotter.plot_total_peptides(st)
        self.assertAlmostEqual(sum(self._bar(fig).y), 100.0, places=3)


class TestBySampleSignificanceAndRelativeNote(unittest.TestCase):
    """Stats must show in By-Sample orientation; relative metric shows a note."""

    GROUPS = ['Ctrl', 'LowDose', 'HighDose']

    def _state(self, **extra):
        rng = np.random.default_rng(6)
        reps = {g: [f'{g}_{i}' for i in (1, 2, 3)] for g in self.GROUPS}
        rows = []
        for i in range(12):
            prot = ['P1', 'P2'][i % 2]
            row = {'Unique Peptide ID': f'pep{i}', 'Protein': prot, 'function': 'ACE-inhibitory'}
            for g in self.GROUPS:
                base = {'Ctrl': 100, 'LowDose': 120, 'HighDose': 400}[g]
                for rc in reps[g]:
                    row[rc] = float(base + rng.normal(0, 8))
                row[f'Avg_{g}'] = float(np.mean([row[rc] for rc in reps[g]]))
            rows.append(row)
        pdct = {'P1': {'name': 'Prot-1'}, 'P2': {'name': 'Prot-2'}}
        params = dict(selected_groups=self.GROUPS, plot_type='Grouped Bar Plots',
                      selected_proteins=['P1', 'P2'],
                      selected_functions=['All Functional Peptides'],
                      abs_or_count='Abundance', metric_type='Absolute',
                      plot_filter='Selected Protein(s)', orientation='By Sample',
                      show_significance=True, log_transform=False)
        params.update(extra)
        st = DataAnalysisState(pd.DataFrame(rows), reps, pdct, params)
        st.run_pipeline()
        return st

    @staticmethod
    def _caption(fig):
        try:
            t = fig.layout.xaxis.title.text or ''
        except Exception:
            t = ''
        return t if 'color:#555555' in t else None

    def test_bysample_selected_protein_shows_stats(self):
        # The reported bug: By Protein filter + By Sample orientation drew no stats.
        st = self._state(orientation='By Sample')
        fig = plotter.create_grouped_bar_plot(st)
        # A significance comparison ran -> the honesty caption naming the test is present.
        cap = self._caption(fig)
        self.assertIsNotNone(cap)
        self.assertIn('Tukey', cap)
        # 3 groups -> a compact-letter scatter-text trace was added.
        text_traces = [t for t in fig.data if getattr(t, 'mode', None) == 'text']
        self.assertGreaterEqual(len(text_traces), 1)

    def test_bysample_matches_byprotein_availability(self):
        # Stats availability must not depend on orientation (the user's point).
        for orient in ('By Sample', 'By Protein'):
            st = self._state(orientation=orient)
            fig = plotter.create_grouped_bar_plot(st)
            self.assertIsNotNone(self._caption(fig), f'{orient}: expected stats caption')

    def test_relative_metric_shows_explanatory_note(self):
        st = self._state(orientation='By Sample', metric_type='Relative')
        fig = plotter.create_grouped_bar_plot(st)
        cap = self._caption(fig)
        self.assertIsNotNone(cap)
        self.assertIn('relative', cap.lower())
        self.assertIn('Absolute', cap)

    def test_count_metric_is_supported(self):
        st = self._state(orientation='By Sample', abs_or_count='Count')
        fig = plotter.create_grouped_bar_plot(st)
        self.assertIsNotNone(self._caption(fig))


class TestDynamicTitles(unittest.TestCase):
    """Tasks 3 & 4: metric-qualified distribution titles; correlation title."""

    def _state(self, **params):
        df = pd.DataFrame({'Unique Peptide ID': ['p'], 'Protein': ['P1'],
                           'Avg_A': [1.0], 'Avg_B': [2.0]})
        base = dict(selected_groups=['A', 'B'], abs_or_count='Abundance',
                    metric_type='Absolute', orientation='By Sample')
        base.update(params)
        return DataAnalysisState(df, {'A': ['A_1'], 'B': ['B_1']}, {}, base)

    def test_absolute_abundance_by_sample_drops_qualifier(self):
        st = self._state(abs_or_count='Abundance', metric_type='Absolute', orientation='By Sample')
        self.assertEqual(plotter._make_title(st), 'Absolute Abundance Distribution')

    def test_relative_abundance_by_function(self):
        st = self._state(abs_or_count='Abundance', metric_type='Relative', orientation='By Function')
        self.assertEqual(plotter._make_title(st), 'Relative Abundance Distribution by Function')

    def test_absolute_count_by_protein(self):
        st = self._state(abs_or_count='Count', metric_type='Absolute', orientation='By Protein')
        self.assertEqual(plotter._make_title(st), 'Absolute Count Distribution by Protein')

    def test_relative_count_by_sample(self):
        st = self._state(abs_or_count='Count', metric_type='Relative', orientation='By Sample')
        self.assertEqual(plotter._make_title(st), 'Relative Count Distribution')

    def test_correlation_title(self):
        st = self._state(plot_type='Corr. Scatter Plots')
        self.assertEqual(plotter._make_title(st, kind='correlation'),
                         'Pairwise Abundance Correlation')

    def test_user_title_overrides(self):
        st = self._state(plot_title='My Custom Title')
        self.assertEqual(plotter._make_title(st), 'My Custom Title')
        self.assertEqual(plotter._make_title(st, kind='correlation'), 'My Custom Title')


class TestCohensD(unittest.TestCase):
    """Two-group effect sizes for the heatmap differential-comparison track."""

    def test_sign_positive_when_a_higher(self):
        # a clearly higher than b -> positive d
        d = stats.cohens_d([10, 11, 12], [1, 2, 3])
        self.assertGreater(d, 0)

    def test_sign_negative_when_b_higher(self):
        d = stats.cohens_d([1, 2, 3], [10, 11, 12])
        self.assertLess(d, 0)

    def test_known_value(self):
        # a=[2,4,6] (mean 4, var 4), b=[1,2,3] (mean 2, var 1)
        # s_pooled = sqrt((2*4 + 2*1)/4) = sqrt(2.5); d = 2/sqrt(2.5)
        d = stats.cohens_d([2, 4, 6], [1, 2, 3])
        self.assertAlmostEqual(d, 2.0 / np.sqrt(2.5), places=10)

    def test_drops_nan(self):
        # NaN dropped listwise; result identical to the clean vectors
        d_nan = stats.cohens_d([2, 4, 6, np.nan], [1, 2, np.nan, 3])
        d_clean = stats.cohens_d([2, 4, 6], [1, 2, 3])
        self.assertAlmostEqual(d_nan, d_clean, places=12)

    def test_undefined_when_too_few(self):
        self.assertTrue(np.isnan(stats.cohens_d([5], [1, 2, 3])))
        self.assertTrue(np.isnan(stats.cohens_d([np.nan, np.nan], [1, 2, 3])))

    def test_undefined_when_zero_pooled_sd(self):
        # both groups internally constant -> pooled SD 0 -> nan (not +inf)
        self.assertTrue(np.isnan(stats.cohens_d([5, 5, 5], [3, 3, 3])))


class TestLog2FoldChange(unittest.TestCase):
    def test_known_value(self):
        # means 8 and 2 -> log2(4) = 2
        self.assertAlmostEqual(stats.log2_fold_change([8, 8], [2, 2]), 2.0, places=12)

    def test_sign(self):
        self.assertLess(stats.log2_fold_change([2, 2], [8, 8]), 0)

    def test_zero_mean_without_eps_is_nan(self):
        self.assertTrue(np.isnan(stats.log2_fold_change([0, 0], [4, 4])))

    def test_eps_pseudocount_defines_ratio(self):
        self.assertTrue(np.isfinite(stats.log2_fold_change([0, 0], [4, 4], eps=1.0)))

    def test_drops_nan(self):
        self.assertAlmostEqual(
            stats.log2_fold_change([8, 8, np.nan], [2, np.nan, 2]), 2.0, places=12)


class TestEffectSizeLabel(unittest.TestCase):
    def test_thresholds(self):
        self.assertEqual(stats.effect_size_label(0.10), 'very small')
        self.assertEqual(stats.effect_size_label(0.30), 'small')
        self.assertEqual(stats.effect_size_label(0.54), 'medium')   # Cheddar β 60-68
        self.assertEqual(stats.effect_size_label(0.90), 'large')
        self.assertEqual(stats.effect_size_label(1.50), 'very large')
        self.assertEqual(stats.effect_size_label(2.50), 'huge')

    def test_uses_magnitude(self):
        self.assertEqual(stats.effect_size_label(-2.5), 'huge')

    def test_nan(self):
        self.assertEqual(stats.effect_size_label(float('nan')), 'n/a')
        self.assertEqual(stats.effect_size_label(None), 'n/a')


class TestSelectorReplicateFlags(unittest.TestCase):
    """get_selector_options exposes replicate presence so the Data Analysis
    dashboard can show the same 'replicate data detected / disabled' banner as
    the Heatmap dashboard (R2-2d Step 2)."""

    def _opts(self, df):
        from data_analysis.services import data_processor as dp
        gdd, renamed, _warn = dp.process_group_data(df)
        pdict = dp.extract_protein_dict(renamed)
        return dp.get_selector_options(renamed, gdd, pdict)

    def test_grouped_columns_flagged(self):
        df = pd.DataFrame({
            'Protein': ['P1', 'P1'], 'Unique Peptide ID': ['a', 'b'],
            "S1 'Grouped: (Bitter)'": [10, 5], "S2 'Grouped: (Bitter)'": [12, 6],
            "S3 'Grouped: (Plain)'": [3, 8],  # single replicate -> not "replicate-level"
            'Avg_Bitter': [11.0, 5.5], 'Avg_Plain': [3.0, 8.0],
        })
        opts = self._opts(df)
        self.assertTrue(opts['has_replicates'])
        self.assertTrue(opts['var_replicates']['Bitter'])       # 2 replicates
        self.assertFalse(opts['var_replicates']['Plain'])       # only 1 replicate

    def test_single_average_file_has_no_replicates(self):
        df = pd.DataFrame({
            'Protein': ['P1'], 'Unique Peptide ID': ['a'],
            'Avg_Bitter': [11.0], 'Avg_Plain': [3.0],
        })
        opts = self._opts(df)
        self.assertFalse(opts['has_replicates'])
        self.assertFalse(any(opts['var_replicates'].values()))


class TestStackedBarTickLabels(unittest.TestCase):
    """Notebook parity: on a stacked bar, an absolute/scaled segment's y-value is
    a proportional share of a group total, not a directly-readable number — so
    y-axis tick labels are hidden unless the metric is Relative (%). A plain,
    non-stacked summary bar (plot_total_peptides) has no such caveat and always
    keeps its tick labels. See data_analysis.ipynb Plotter.plot_stacked_bar_scaled
    (`showticklabels_tf`)."""

    GROUPS = ['Ctrl', 'LowDose', 'HighDose']

    def _dataset(self):
        rng = np.random.default_rng(0)
        reps = {g: [f'{g}_{i}' for i in range(1, 4)] for g in self.GROUPS}
        rows = []
        for i in range(12):
            fn = 'ACE-inhibitory' if i % 2 == 0 else 'Antioxidant'
            prot = 'P02666' if i < 6 else 'P02662'
            row = {'Unique Peptide ID': f'pep{i}', 'Protein': prot, 'function': fn}
            for g in self.GROUPS:
                base = {'Ctrl': 100, 'LowDose': 120, 'HighDose': 500}[g]
                for rc in reps[g]:
                    row[rc] = float(base + rng.normal(0, 8))
                row[f'Avg_{g}'] = float(np.mean([row[rc] for rc in reps[g]]))
            rows.append(row)
        protein_dict = {'P02666': {'name': 'Beta-casein'},
                        'P02662': {'name': 'Alpha-S1-casein'}}
        return pd.DataFrame(rows), reps, protein_dict

    def _state(self, reps, merged, protein_dict, **extra):
        params = dict(selected_groups=self.GROUPS,
                      selected_functions=['All Functional Peptides'],
                      selected_proteins=['All Proteins (No Filter)'],
                      abs_or_count='Abundance', metric_type='Absolute')
        params.update(extra)
        st = DataAnalysisState(merged, reps, protein_dict, params)
        st.run_pipeline()
        return st

    def _showticklabels(self, fig):
        return fig.layout.yaxis.showticklabels

    def test_no_filter_stacked_hides_ticks_for_absolute_shows_for_relative(self):
        merged, reps, pdct = self._dataset()
        st_abs = self._state(reps, merged, pdct, orientation='By Sample',
                              plot_filter='No Filter', metric_type='Absolute')
        fig_abs = plotter.plot_stacked_bar_scaled(st_abs)
        self.assertIs(self._showticklabels(fig_abs), False,
                       'scaled/absolute stacked-by-sample bar should hide y ticks')

        st_rel = self._state(reps, merged, pdct, orientation='By Sample',
                              plot_filter='No Filter', metric_type='Relative')
        fig_rel = plotter.plot_stacked_bar_scaled(st_rel)
        self.assertTrue(self._showticklabels(fig_rel),
                        'relative (%) stacked-by-sample bar should show y ticks')

    def test_selected_function_stacked_hides_ticks_for_absolute_shows_for_relative(self):
        merged, reps, pdct = self._dataset()
        st_abs = self._state(reps, merged, pdct, orientation='By Sample',
                              plot_filter='Selected Function(s)', metric_type='Absolute',
                              selected_functions=['ACE-inhibitory', 'Antioxidant'])
        fig_abs = plotter.plot_stacked_bar_scaled(st_abs)
        self.assertIs(self._showticklabels(fig_abs), False,
                       'scaled/absolute stacked-by-function bar should hide y ticks')

        st_rel = self._state(reps, merged, pdct, orientation='By Sample',
                              plot_filter='Selected Function(s)', metric_type='Relative',
                              selected_functions=['ACE-inhibitory', 'Antioxidant'])
        fig_rel = plotter.plot_stacked_bar_scaled(st_rel)
        self.assertTrue(self._showticklabels(fig_rel),
                        'relative (%) stacked-by-function bar should show y ticks')

    def test_plain_summary_plot_always_shows_ticks(self):
        # plot_total_peptides: one non-stacked bar per sample, no scaling, no
        # layered segments -> tick labels are always meaningful and shown.
        merged, reps, pdct = self._dataset()
        st = self._state(reps, merged, pdct, orientation='By Sample',
                         plot_filter='No Filter', metric_type='Absolute')
        fig = plotter.plot_total_peptides(st)
        self.assertIsNot(self._showticklabels(fig), False,
                         'plain summary bar should not hide y ticks')


if __name__ == '__main__':
    unittest.main()
