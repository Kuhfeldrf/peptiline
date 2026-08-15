"""
Tests for the heatmap differential-comparison feature
(manuscript/docs/HEATMAP_DIFFERENTIAL_STATS.md).

Step 1 — replicate-aware loader: the heatmap data_processor must parse the
``"<col> 'Grouped: (Group)'"`` replicate framework and surface per-group
replicate columns so a pooled SD (Cohen's d) can be computed downstream.

Run with:
    cd /home/kuhfeldrf/mbpdb/include/peptide
    /home/kuhfeldrf/mbpdb/.venv/bin/python -m pytest tests/test_heatmap_differential.py -v
"""
import io
import json
import os
import unittest

import numpy as np
import pandas as pd

# conftest.py bootstraps Django; import the services directly.
from heatmap_viz.services import data_processor as hdp
from heatmap_viz.services import heatmap_renderer as hdp_render


def _csv_bytes(df: pd.DataFrame) -> io.BytesIO:
    buf = io.BytesIO()
    df.to_csv(buf, index=False)
    buf.seek(0)
    return buf


class TestParseGroupedReplicates(unittest.TestCase):
    def test_parses_and_renames(self):
        df = pd.DataFrame({
            'Protein': ['P02666'],
            "S1 'Grouped: (Bitter)'": [10],
            "S2 'Grouped: (Bitter)'": [12],
            "S3 'Grouped: (NonBitter)'": [3],
        })
        renamed, reps = hdp._parse_grouped_replicates(df)
        self.assertEqual(reps, {'Bitter': ['S1', 'S2'], 'NonBitter': ['S3']})
        # Columns renamed back to base names for indexing by replicate.
        self.assertIn('S1', renamed.columns)
        self.assertNotIn("S1 'Grouped: (Bitter)'", renamed.columns)

    def test_multi_group_membership(self):
        # A replicate can belong to more than one group (semicolon list).
        df = pd.DataFrame({"S1 'Grouped: (A; B)'": [1]})
        _, reps = hdp._parse_grouped_replicates(df)
        self.assertEqual(reps, {'A': ['S1'], 'B': ['S1']})

    def test_no_grouped_columns_is_noop(self):
        df = pd.DataFrame({'Protein': ['P1'], 'Avg_X': [1.0]})
        renamed, reps = hdp._parse_grouped_replicates(df)
        self.assertEqual(reps, {})
        self.assertTrue(renamed.equals(df))


class TestLoadMergedFileReplicates(unittest.TestCase):
    def _df(self):
        return pd.DataFrame({
            'Protein': ['P02666', 'P02666'],
            'start': [60, 73],
            'end': [68, 79],
            'Unique Peptide ID': ['pep1', 'pep2'],
            "S1 'Grouped: (Bitter)'": [10, 5],
            "S2 'Grouped: (Bitter)'": [12, 6],
            "S3 'Grouped: (Bitter)'": [11, 4],
            "S4 'Grouped: (NonBitter)'": [3, 8],
            "S5 'Grouped: (NonBitter)'": [2, 9],
            "S6 'Grouped: (NonBitter)'": [4, 7],
            'Avg_Bitter': [11.0, 5.0],
            'Avg_NonBitter': [3.0, 8.0],
        })

    def test_replicate_columns_attached(self):
        df, gdd, pdict, col_order, err = hdp.load_merged_file(
            _csv_bytes(self._df()), 'merged.csv')
        self.assertIsNone(err)
        # group_data_dict keyed by ordinal; Bitter first, NonBitter second.
        by_var = {v['grouping_variable']: v for v in gdd.values()}
        self.assertEqual(by_var['Bitter']['replicate_columns'], ['S1', 'S2', 'S3'])
        self.assertEqual(by_var['NonBitter']['replicate_columns'], ['S4', 'S5', 'S6'])
        self.assertEqual(sorted(col_order), ['Bitter', 'NonBitter'])

    def test_single_average_file_disables_comparison(self):
        # No 'Grouped:' columns -> replicate_columns empty, load still succeeds.
        df = pd.DataFrame({
            'Protein': ['P02666'], 'start': [60], 'end': [68],
            'Unique Peptide ID': ['pep1'],
            'Avg_Bitter': [11.0], 'Avg_NonBitter': [3.0],
        })
        _, gdd, _, _, err = hdp.load_merged_file(_csv_bytes(df), 'merged.csv')
        self.assertIsNone(err)
        for v in gdd.values():
            self.assertEqual(v['replicate_columns'], [])


class TestSelectorOptionsReplicates(unittest.TestCase):
    """Step 2/6 — get_selector_options exposes which variables are replicate-
    backed so the UI can gate the comparison controls + banner."""

    def test_exposes_replicate_flags(self):
        df = pd.DataFrame({
            'Protein': ['P02666', 'P02666'],
            'start': [60, 73], 'end': [68, 79],
            "S1 'Grouped: (Bitter)'": [10, 5],
            "S2 'Grouped: (Bitter)'": [12, 6],
            'Avg_Bitter': [11.0, 5.0],
            'Avg_Plain': [3.0, 8.0],  # no replicate columns
        })
        merged, gdd, pdict, col_order, err = hdp.load_merged_file(
            _csv_bytes(df), 'merged.csv')
        self.assertIsNone(err)
        opts = hdp.get_selector_options(merged, gdd, pdict, col_order)
        self.assertTrue(opts['has_replicates'])
        self.assertTrue(opts['var_replicates']['Bitter'])
        self.assertFalse(opts['var_replicates']['Plain'])

    def test_single_average_file_no_replicates(self):
        df = pd.DataFrame({
            'Protein': ['P02666'], 'start': [60], 'end': [68],
            'Avg_Bitter': [11.0], 'Avg_Plain': [3.0],
        })
        merged, gdd, pdict, col_order, err = hdp.load_merged_file(
            _csv_bytes(df), 'merged.csv')
        opts = hdp.get_selector_options(merged, gdd, pdict, col_order)
        self.assertFalse(opts['has_replicates'])
        self.assertFalse(any(opts['var_replicates'].values()))


class TestContrastTrack(unittest.TestCase):
    """Step 4 — per-position contrast builder over replicate position-means."""

    def _peptides(self):
        # Two peptides spanning a 10-residue protein; 3 replicates per group.
        # Bitter reps clearly higher than NonBitter at every covered position.
        return pd.DataFrame({
            'start': [1, 4],
            'end': [5, 8],
            'B1': [10, 20], 'B2': [12, 22], 'B3': [11, 21],
            'N1': [3, 5], 'N2': [2, 6], 'N3': [4, 4],
        })

    def test_position_means_shape_and_drop(self):
        seq = 'ACDEFGHIKL'  # length 10
        m = hdp_render.calculate_replicate_position_means(seq, self._peptides(), ['B1', 'B2', 'B3'])
        self.assertEqual(list(m.columns), ['B1', 'B2', 'B3'])
        self.assertEqual(len(m), 10)
        # position 1 (index 0) covered only by peptide 1 -> B1 == 10
        self.assertAlmostEqual(m.iloc[0]['B1'], 10.0)
        # positions 9-10 (index 8,9) covered by neither peptide -> NaN (dropped)
        self.assertTrue(m.iloc[8].isna().all())

    def test_smd_positive_when_a_higher(self):
        seq = 'ACDEFGHIKL'
        peps = self._peptides()
        a = hdp_render.calculate_replicate_position_means(seq, peps, ['B1', 'B2', 'B3'])
        b = hdp_render.calculate_replicate_position_means(seq, peps, ['N1', 'N2', 'N3'])
        track = hdp_render.calculate_contrast_track(a, b, metric='smd')
        self.assertEqual(len(track), 10)
        # covered positions: Bitter >> NonBitter -> strongly positive d
        self.assertGreater(track[0], 1.0)
        # uncovered tail -> NaN gap, not a value
        self.assertTrue(np.isnan(track[8]))

    def test_log2fc_metric(self):
        seq = 'ACDEFGHIKL'
        peps = self._peptides()
        a = hdp_render.calculate_replicate_position_means(seq, peps, ['B1', 'B2', 'B3'])
        b = hdp_render.calculate_replicate_position_means(seq, peps, ['N1', 'N2', 'N3'])
        track = hdp_render.calculate_contrast_track(a, b, metric='log2fc')
        self.assertGreater(track[0], 0)  # A higher -> positive log2FC

    def test_unequal_length_alignment(self):
        # Shorter B frame -> positions beyond it are NaN, no crash.
        a = pd.DataFrame({'r1': [1.0, 2.0, 3.0], 'r2': [1.1, 2.1, 3.1], 'r3': [0.9, 1.9, 2.9]})
        b = pd.DataFrame({'s1': [0.1, 0.2], 's2': [0.1, 0.2], 's3': [0.1, 0.2]})
        track = hdp_render.calculate_contrast_track(a, b, metric='smd')
        self.assertEqual(len(track), 3)
        self.assertTrue(np.isnan(track[2]))


class TestSignedYTicks(unittest.TestCase):
    """Step 5 — symmetric diverging axis for the SMD / log2FC track."""

    def test_symmetric_about_zero(self):
        ticks = hdp_render.calculate_signed_y_ticks(-0.3, 0.54)
        self.assertEqual(ticks[1], 0.0)
        self.assertEqual(ticks[0], -ticks[2])

    def test_bound_covers_max_magnitude(self):
        # bound must be >= the largest |value| seen
        for lo, hi in [(-0.3, 0.54), (-3.49, 1.2), (-0.1, 0.05)]:
            ticks = hdp_render.calculate_signed_y_ticks(lo, hi)
            self.assertGreaterEqual(ticks[2], max(abs(lo), abs(hi)))

    def test_negative_dominant(self):
        ticks = hdp_render.calculate_signed_y_ticks(-3.49, 0.2)
        self.assertGreaterEqual(ticks[2], 3.49)
        self.assertEqual(ticks[0], -ticks[2])

    def test_degenerate_and_nan(self):
        self.assertEqual(hdp_render.calculate_signed_y_ticks(0, 0), [-1.0, 0.0, 1.0])
        self.assertEqual(
            hdp_render.calculate_signed_y_ticks(float('nan'), float('nan')),
            [-1.0, 0.0, 1.0])


class TestRenderBranchComparison(unittest.TestCase):
    """Step 5 — render branch: update_plot wires the contrast track into the
    interactive figure (signed line replacing the abundance line) when
    comparison_mode is on, and validates the series otherwise."""

    SEQ = 'MKVLILACDEFGHIKLMNPQRSTVWY'  # 26-residue toy protein

    def _available(self):
        """Run the real loader → build_available_data_variables pipeline so the
        variable dict carries replicate_columns + peptide_df (as in production)."""
        df = pd.DataFrame({
            'Protein': ['P02666', 'P02666'],
            'protein_name': ['Beta-casein', 'Beta-casein'],
            'protein_species': ['Bos taurus', 'Bos taurus'],
            'start': [3, 10],
            'end': [9, 16],
            'function': ['bitter', None],
            'Unique Peptide ID': ['pep1', 'pep2'],
            "S1 'Grouped: (Bitter)'": [10, 20],
            "S2 'Grouped: (Bitter)'": [12, 22],
            "S3 'Grouped: (Bitter)'": [11, 21],
            "S4 'Grouped: (NonBitter)'": [3, 5],
            "S5 'Grouped: (NonBitter)'": [2, 6],
            "S6 'Grouped: (NonBitter)'": [4, 4],
            'Avg_Bitter': [11.0, 21.0],
            'Avg_NonBitter': [3.0, 5.0],
        })
        merged, gdd, pdict, _col_order, err = hdp.load_merged_file(
            _csv_bytes(df), 'merged.csv')
        self.assertIsNone(err)
        pdict['P02666']['sequence'] = self.SEQ  # inject sequence (no FASTA in test)
        avail, _msgs = hdp.build_available_data_variables(
            merged, pdict, gdd, ['P02666'], ['Bitter', 'NonBitter'])
        return avail

    def _update(self, avail, **overrides):
        params = dict(
            available_data_variables_dict=avail,
            ms_average_choice='yes', bio_or_pep='no',
            selected_peptides=[], selected_functions=[],
            hm_selected_color='RdYlGn_r', lp_selected_color='Set3',
            avglp_selected_color='Dark2',
            xaxis_label='', yaxis_label='', yaxis_position=5,
            legend_title_input_1='Sample Type:', legend_title_input_2='Peptide Counts:',
            legend_title_input_3='Abundance:',
            filter_type='all-peptides', log_transform=False,
            manual_y_axis=False, y_min_manual=0.0, y_max_manual=1.0,
            plot_landscape_interactive=True,
        )
        params.update(overrides)
        return hdp_render.update_plot(**params)

    def test_contrast_line_replaces_abundance(self):
        avail = self._available()
        _p, _l, _c, _fig_port, fig_land, errors, _notes = self._update(
            avail, comparison_mode=True,
            series_a='P02666_Bitter', series_b='P02666_NonBitter',
            comparison_metric='smd')
        self.assertFalse([e for e in errors if 'Comparison' in e], errors)
        self.assertIsNotNone(fig_land)
        # Line-panel traces are Scatter(mode='lines'); density tile rows are Bar
        # (still named after each group — that is the heatmap the contrast sits on).
        line_names = [
            t.name or '' for t in fig_land.data
            if t.type == 'scatter' and 'lines' in (t.mode or '')
        ]
        # the only line trace is the contrast ("A vs B"); no abundance line per group
        self.assertTrue(any(' vs ' in n for n in line_names), line_names)
        self.assertFalse(any(n in ('Bitter', 'NonBitter') for n in line_names), line_names)
        # its legend group header is "Comparison", not the abundance "Sample Type:"
        contrast = next(t for t in fig_land.data
                        if t.type == 'scatter' and 'lines' in (t.mode or ''))
        self.assertEqual(contrast.legendgrouptitle.text, 'Comparison')

    def test_signed_axis_symmetric_with_zeroline(self):
        avail = self._available()
        *_rest, fig_land, _errors, _notes = self._update(
            avail, comparison_mode=True,
            series_a='P02666_Bitter', series_b='P02666_NonBitter',
            comparison_metric='smd')
        lay = fig_land.layout.to_plotly_json()
        # the line-panel y-axis (yaxis, first row) carries the zero baseline and a
        # range symmetric about 0.
        yax = lay.get('yaxis', {})
        self.assertTrue(yax.get('zeroline'))
        rng = yax.get('range')
        self.assertIsNotNone(rng)
        self.assertAlmostEqual(rng[0], -rng[1], places=6)

    def test_log2fc_metric_branch(self):
        avail = self._available()
        *_rest, fig_land, errors, _notes = self._update(
            avail, comparison_mode=True,
            series_a='P02666_Bitter', series_b='P02666_NonBitter',
            comparison_metric='log2fc')
        self.assertFalse([e for e in errors if 'Comparison' in e], errors)
        # metric name lives on the y-axis annotation now, not the trace name
        ann = [a.get('text', '') for a in fig_land.layout.to_plotly_json().get('annotations', [])]
        self.assertTrue(any('Fold Change' in t for t in ann), ann)

    def test_smd_axis_title_spelled_out(self):
        avail = self._available()
        *_rest, fig_land, _errors, _notes = self._update(
            avail, comparison_mode=True,
            series_a='P02666_Bitter', series_b='P02666_NonBitter',
            comparison_metric='smd')
        ann = [a.get('text', '') for a in fig_land.layout.to_plotly_json().get('annotations', [])]
        self.assertTrue(any('Standardized Mean Difference' in t for t in ann), ann)

    def test_same_series_rejected(self):
        avail = self._available()
        *_rest, errors, _notes = self._update(
            avail, comparison_mode=True,
            series_a='P02666_Bitter', series_b='P02666_Bitter',
            comparison_metric='smd')
        self.assertTrue(any('distinct' in e for e in errors), errors)

    def test_comparison_off_keeps_abundance_line(self):
        avail = self._available()
        *_rest, fig_land, _errors, _notes = self._update(avail, comparison_mode=False)
        names = [t.name or '' for t in fig_land.data]
        # abundance lines present, no contrast trace
        self.assertFalse(any('SMD' in n for n in names), names)
        self.assertTrue(any(n in ('Bitter', 'NonBitter') for n in names), names)

    def _available_three_vars(self):
        """One protein, THREE replicate-backed variables — to prove the density
        heatmap is restricted to only the two compared series."""
        df = pd.DataFrame({
            'Protein': ['P02666', 'P02666'],
            'protein_name': ['Beta-casein', 'Beta-casein'],
            'protein_species': ['Bos taurus', 'Bos taurus'],
            'start': [3, 10], 'end': [9, 16],
            'function': ['bitter', None],
            'Unique Peptide ID': ['pep1', 'pep2'],
            "B1 'Grouped: (Bitter)'": [10, 20], "B2 'Grouped: (Bitter)'": [12, 22],
            "N1 'Grouped: (NonBitter)'": [3, 5], "N2 'Grouped: (NonBitter)'": [2, 6],
            "P1 'Grouped: (Plain)'": [7, 8], "P2 'Grouped: (Plain)'": [6, 9],
            'Avg_Bitter': [11.0, 21.0], 'Avg_NonBitter': [3.0, 5.0], 'Avg_Plain': [6.5, 8.5],
        })
        merged, gdd, pdict, _c, err = hdp.load_merged_file(_csv_bytes(df), 'm.csv')
        self.assertIsNone(err)
        pdict['P02666']['sequence'] = self.SEQ
        avail, _m = hdp.build_available_data_variables(
            merged, pdict, gdd, ['P02666'], ['Bitter', 'NonBitter', 'Plain'])
        return avail

    def test_heatmap_restricted_to_compared_series(self):
        avail = self._available_three_vars()
        self.assertEqual(len(avail), 3)  # three variables selected
        *_rest, fig_land, errors, _notes = self._update(
            avail, comparison_mode=True,
            series_a='P02666_Bitter', series_b='P02666_NonBitter',
            comparison_metric='smd')
        self.assertFalse([e for e in errors if 'Comparison' in e], errors)
        # density tile rows are Bar traces named after the variable; only the two
        # compared series (Bitter, NonBitter) should appear — NOT the third (Plain).
        bar_names = {t.name for t in fig_land.data if t.type == 'bar'}
        self.assertIn('Bitter', bar_names)
        self.assertIn('NonBitter', bar_names)
        self.assertNotIn('Plain', bar_names)

    def _available_two_proteins(self):
        """Two proteins, same variable — the A1/A2 style comparison. The density
        rows must be protein-labeled so they are distinguishable."""
        df = pd.DataFrame({
            'Protein': ['P02666', 'P02662'],
            'protein_name': ['Beta-casein', 'Alpha-S1-casein'],
            'protein_species': ['Bos taurus', 'Bos taurus'],
            'start': [3, 3], 'end': [9, 9],
            'function': ['bitter', 'bitter'],
            'Unique Peptide ID': ['pep1', 'pep2'],
            "B1 'Grouped: (Bitter)'": [10, 4], "B2 'Grouped: (Bitter)'": [12, 5],
            'Avg_Bitter': [11.0, 4.5],
        })
        merged, gdd, pdict, _c, err = hdp.load_merged_file(_csv_bytes(df), 'm.csv')
        self.assertIsNone(err)
        pdict['P02666']['sequence'] = self.SEQ
        pdict['P02662']['sequence'] = self.SEQ
        avail, _m = hdp.build_available_data_variables(
            merged, pdict, gdd, ['P02666', 'P02662'], ['Bitter'])
        return avail

    def test_two_protein_labels_include_protein_name(self):
        avail = self._available_two_proteins()
        *_rest, fig_land, errors, _notes = self._update(
            avail, comparison_mode=True,
            series_a='P02666_Bitter', series_b='P02662_Bitter',
            comparison_metric='smd')
        self.assertFalse([e for e in errors if 'Comparison' in e], errors)
        # contrast trace name distinguishes the two proteins (not "Bitter vs Bitter")
        line_names = [t.name or '' for t in fig_land.data
                      if t.type == 'scatter' and 'lines' in (t.mode or '')]
        self.assertTrue(line_names)
        self.assertIn('Beta-casein', line_names[0])
        self.assertIn('Alpha-S1-casein', line_names[0])
        # density rows are protein-labeled and distinct
        bar_names = {t.name for t in fig_land.data if t.type == 'bar'}
        self.assertIn('Beta-casein', bar_names)
        self.assertIn('Alpha-S1-casein', bar_names)

    def test_multi_protein_label_prefixes_protein_name(self):
        """Comparing proteins: the bare sample name ('Bitter') repeats across
        proteins, so the display label must carry the protein name while
        `var_label` keeps the bare sample name for downstream callers."""
        avail = self._available_two_proteins()
        self.assertEqual(
            {vd['label'] for vd in avail.values()},
            {'Beta-casein Bitter', 'Alpha-S1-casein Bitter'})
        self.assertEqual({vd['var_label'] for vd in avail.values()}, {'Bitter'})

    def test_multi_protein_labels_in_row_and_legend(self):
        """Without comparison mode, each protein's density row + averaged-line
        legend entry is protein-labeled so the two proteins are distinguishable."""
        avail = self._available_two_proteins()
        *_rest, fig_land, errors, _notes = self._update(avail)
        self.assertFalse(errors, errors)
        bar_names = {t.name for t in fig_land.data if t.type == 'bar'}
        self.assertIn('Beta-casein Bitter', bar_names)
        self.assertIn('Alpha-S1-casein Bitter', bar_names)
        line_names = {t.name for t in fig_land.data
                      if t.type == 'scatter' and 'lines' in (t.mode or '')}
        self.assertIn('Beta-casein Bitter', line_names)
        self.assertIn('Alpha-S1-casein Bitter', line_names)

    def test_single_protein_keeps_bare_sample_names(self):
        """Comparing samples on one protein: the sample name alone is
        unambiguous and must stay in the row/legend label (no protein prefix)."""
        avail = self._available_three_vars()  # single protein, three samples
        self.assertEqual(
            {vd['label'] for vd in avail.values()},
            {'Bitter', 'NonBitter', 'Plain'})


class TestTransferFromDtFasta(unittest.TestCase):
    """Transfer from Data Transformation should carry an optional FASTA into the
    heatmap so protein sequences are loaded by default, and expose its filename
    so the page can display it (mirrors a manual FASTA upload)."""

    def _df(self):
        return pd.DataFrame({
            'Protein': ['P02666', 'P02666'],
            'start': [60, 73], 'end': [68, 79],
            'Unique Peptide ID': ['pep1', 'pep2'],
            'Avg_Bitter': [11.0, 5.0], 'Avg_NonBitter': [3.0, 8.0],
        })

    def setUp(self):
        import tempfile
        from django.test import Client
        self.client = Client()
        self.dt_dir = tempfile.mkdtemp(prefix='dt_test_')
        self.hm_dirs = []

    def tearDown(self):
        import shutil
        shutil.rmtree(self.dt_dir, ignore_errors=True)
        for d in self.hm_dirs:
            shutil.rmtree(d, ignore_errors=True)

    def _seed_dt_session(self, with_fasta):
        import json as _json
        self._df().to_pickle(os.path.join(self.dt_dir, 'merged_df.pkl'))
        if with_fasta:
            with open(os.path.join(self.dt_dir, 'fasta_file.orig'), 'w') as fh:
                fh.write('>P02666 Beta-casein\nRELEELNVPGEIVESLSSSEESITR\n')
            with open(os.path.join(self.dt_dir, 'uploaded_file_names.json'), 'w') as fh:
                _json.dump({'fasta_file': 'my_seqs.fasta'}, fh)
        session = self.client.session
        session['dt_work_dir'] = self.dt_dir
        session.save()

    def _post_transfer(self):
        resp = self.client.post('/heatmap/transfer-from-dt/')
        hm_dir = self.client.session.get('hm_work_dir')
        if hm_dir:
            self.hm_dirs.append(hm_dir)
        return resp, hm_dir

    def test_fasta_carried_and_named(self):
        self._seed_dt_session(with_fasta=True)
        resp, hm_dir = self._post_transfer()
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertTrue(data.get('success'))
        self.assertEqual(data.get('fasta_filename'), 'my_seqs.fasta')
        # Sequence must have been merged into the heatmap's protein_dict.
        with open(os.path.join(hm_dir, 'protein_dict.json')) as f:
            pdict = json.load(f)
        self.assertIn('sequence', pdict.get('P02666', {}))
        self.assertTrue(pdict['P02666']['sequence'])

    def test_no_fasta_leaves_it_cleared(self):
        # No fasta_file.orig (e.g. after a Data Transformation "Start Over").
        self._seed_dt_session(with_fasta=False)
        resp, _ = self._post_transfer()
        self.assertEqual(resp.status_code, 200)
        self.assertIsNone(resp.json().get('fasta_filename'))


if __name__ == '__main__':
    unittest.main()
