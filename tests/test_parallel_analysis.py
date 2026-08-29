"""Parallel post-load analysis on merged-file upload.

The Data Analysis and Heatmap upload paths run the independent passes over the
peptide table (protein-info extraction, abundance totals, function counts /
replicate parsing) concurrently. These tests pin the contract: the parallel
path returns exactly what the old sequential calls did, and run_tasks handles
the edge cases (empty / single task, exception propagation).
"""
import time
import unittest

import pandas as pd

from utils.parallel import run_tasks


class TestRunTasks(unittest.TestCase):
    def test_empty_and_single(self):
        self.assertEqual(run_tasks({}), {})
        self.assertEqual(run_tasks({'a': lambda: 1}), {'a': 1})

    def test_results_keyed_correctly(self):
        out = run_tasks({'a': lambda: 'A', 'b': lambda: 'B', 'c': lambda: 'C'})
        self.assertEqual(out, {'a': 'A', 'b': 'B', 'c': 'C'})

    def test_exception_propagates(self):
        def boom():
            raise ValueError('nope')
        with self.assertRaises(ValueError):
            run_tasks({'ok': lambda: 1, 'bad': boom})

    def test_runs_concurrently(self):
        # Three 0.15s GIL-releasing sleeps finish well under the 0.45s a
        # sequential run would take.
        start = time.monotonic()
        run_tasks({k: (lambda: time.sleep(0.15)) for k in ('a', 'b', 'c')})
        self.assertLess(time.monotonic() - start, 0.35)


def _sample_df():
    return pd.DataFrame({
        'Unique Peptide ID': [f'p{i}' for i in range(6)],
        'Protein': ['P1', 'P1', 'P2', 'P2;P3', 'P3', 'P1'],
        'protein_name': ['Alpha', 'Alpha', 'Beta', 'Beta', 'Gamma', 'Alpha'],
        'protein_species': ['Bos'] * 6,
        'function': ['ACE-inhibitory; Antioxidant', 'ACE-inhibitory', None,
                     'Antioxidant', 'ACE-inhibitory', None],
        "S1 'Grouped: (Hi)'": [10, 20, 30, 40, 50, 60],
        "S2 'Grouped: (Hi)'": [11, 21, 31, 41, 51, 61],
        'Avg_Hi': [10.5, 20.5, 30.5, 40.5, 50.5, 60.5],
        'Avg_Lo': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    })


class TestDataAnalysisParallelMatchesSequential(unittest.TestCase):
    def test_analyze_merged_dataframe_matches_piecewise(self):
        from data_analysis.services import data_processor as dp
        df = _sample_df()
        gdd, renamed, _ = dp.process_group_data(df)

        seq_pdict = dp.extract_protein_dict(renamed)
        seq_opts = dp.get_selector_options(renamed, gdd, seq_pdict)

        par_pdict, par_opts = dp.analyze_merged_dataframe(renamed, gdd)

        self.assertEqual(par_pdict, seq_pdict)
        self.assertEqual(par_opts, seq_opts)
        # sanity: the aggregates actually populated
        self.assertTrue(par_opts['protein_ids'])
        self.assertIn('ACE-inhibitory', par_opts['functions'])


class TestHeatmapParallelLoad(unittest.TestCase):
    def test_load_merged_file_builds_protein_dict_and_replicates(self):
        import io
        from heatmap_viz.services import data_processor as hdp
        sdf = _sample_df()
        sdf['start'] = [1, 4, 1, 10, 5, 20]
        sdf['end'] = [8, 12, 9, 18, 14, 28]
        csv = sdf.to_csv(index=False).encode()
        df, gdd, pdict, col_order, err = hdp.load_merged_file(io.BytesIO(csv), 'm.csv')
        self.assertIsNone(err)
        self.assertEqual(set(pdict), {'P1', 'P2', 'P2;P3', 'P3'})
        # 'Grouped:' columns were renamed to their base name
        self.assertIn('S1', df.columns)
        hi = next(v for v in gdd.values() if v['grouping_variable'] == 'Hi')
        self.assertEqual(sorted(hi['replicate_columns']), ['S1', 'S2'])


class TestBlankRowFilter(unittest.TestCase):
    """The all-blank trailing-row drop is column-wise vectorised now (a row-wise
    .apply(axis=1) cost ~4.5s of a 5s load on the 40k-row case study). Behaviour
    must be unchanged: rows that are entirely empty / whitespace go, real rows
    stay, NaN-only rows are handled by the preceding dropna."""

    def test_data_analysis_load_file_drops_blank_rows(self):
        import io
        from data_analysis.services import data_processor as dp
        csv = (
            "Unique Peptide ID,Protein,Avg_A\n"
            "p1,P1,10\n"
            " , , \n"          # whitespace-only row
            ",,\n"              # all-empty row
            "p2,P2,20\n"
        )
        df, err = dp.load_file(io.BytesIO(csv.encode()), 'm.csv')
        self.assertIsNone(err)
        self.assertEqual(list(df['Unique Peptide ID']), ['p1', 'p2'])

    def test_heatmap_load_merged_file_drops_blank_rows(self):
        import io
        from heatmap_viz.services import data_processor as hdp
        csv = (
            "Protein,Peptide,start,end,Avg_A\n"
            "P1,AAA,1,3,10\n"
            " , , , , \n"
            "P2,BBB,4,6,20\n"
        )
        df, gdd, pdict, col_order, err = hdp.load_merged_file(io.BytesIO(csv.encode()), 'm.csv')
        self.assertIsNone(err)
        self.assertEqual(len(df), 2)


if __name__ == '__main__':
    unittest.main()
