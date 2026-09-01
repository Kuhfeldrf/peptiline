"""
Tests for the "arrange & rename" panel wiring in Data Analysis:

  * ``selected_groups`` order is honoured on the figure;
  * ``label_overrides`` rename axis / legend / hover text everywhere the name
    appears, keyed by the ORIGINAL value;
  * a rename is display-only — it must not move or recolour a series.

Run with:
    .venv/bin/python -m pytest tests/test_data_analysis_labels.py -v
"""
import json
import unittest

import numpy as np
import pandas as pd

# conftest.py bootstraps Django; import the services directly.
from data_analysis.services import plotter
from data_analysis.services import data_processor as dp


GROUPS = ['Ctrl', 'LowDose', 'HighDose']


def _dataset():
    rng = np.random.default_rng(0)
    reps = {g: [f'{g}_{i}' for i in range(1, 4)] for g in GROUPS}
    rows = []
    for i in range(10):
        prot = 'P02666' if i < 5 else 'P02662'
        row = {'Unique Peptide ID': f'pep{i}', 'Protein': prot,
               'function': 'ACE-inhibitory' if i % 2 == 0 else 'Antioxidant'}
        for g in GROUPS:
            base = {'Ctrl': 100.0, 'LowDose': 160.0, 'HighDose': 500.0}[g]
            for rc in reps[g]:
                row[rc] = float(base + rng.normal(0, 5))
            row[f'Avg_{g}'] = float(np.mean([row[rc] for rc in reps[g]]))
        rows.append(row)
    protein_dict = {'P02666': {'name': 'Beta-casein'},
                    'P02662': {'name': 'Alpha-S1-casein'}}
    return pd.DataFrame(rows), reps, protein_dict


def _plot(**extra):
    merged, gdd, pdict = _dataset()
    params = dict(
        selected_groups=GROUPS,
        selected_proteins=['No Protein Filter'],
        selected_functions=['No Function Filter'],
        plot_type='Grouped Bar Plots', orientation='By Sample',
        abs_or_count='Abundance', metric_type='Absolute',
    )
    params.update(extra)
    fig_json, warnings = plotter.generate_plot(merged, gdd, pdict, params)
    assert fig_json, warnings
    return json.loads(fig_json), warnings


class TestSampleOrder(unittest.TestCase):
    def test_group_order_follows_selection(self):
        fig, _ = _plot(selected_groups=['HighDose', 'Ctrl', 'LowDose'])
        xs = [list(t['x']) for t in fig['data'] if t.get('x')]
        self.assertTrue(xs)
        self.assertEqual(xs[0], ['HighDose', 'Ctrl', 'LowDose'])


class TestRename(unittest.TestCase):
    def test_rename_appears_on_axis_and_hover_not_original(self):
        fig, _ = _plot(
            selected_groups=['Ctrl', 'LowDose', 'HighDose'],
            label_overrides={'sample_groups': {'HighDose': 'High dose (10 mg)'}},
        )
        blob = json.dumps(fig)
        self.assertIn('High dose (10 mg)', blob)

        xs = [list(t['x']) for t in fig['data'] if t.get('x')][0]
        self.assertIn('High dose (10 mg)', xs)
        self.assertNotIn('HighDose', xs)

        # hover text is rewritten too
        hovers = []
        for t in fig['data']:
            ht = t.get('hovertext')
            if isinstance(ht, list):
                hovers.extend(ht)
            elif isinstance(ht, str):
                hovers.append(ht)
        joined = ' '.join(hovers)
        if joined:  # totals plot carries per-bar hovertext
            self.assertIn('High dose (10 mg)', joined)
            self.assertNotIn('HighDose', joined)

    def test_protein_rename_by_accession(self):
        fig, _ = _plot(
            selected_proteins=['P02666', 'P02662'],
            orientation='By Protein', plot_type='Stacked Bar Plots',
            label_overrides={'proteins': {'P02666': 'Casein beta chain'}},
        )
        blob = json.dumps(fig, ensure_ascii=False)
        self.assertIn('Casein beta chain', blob)
        # the accession's resolved name is remapped everywhere it was printed
        xs = [list(t['x']) for t in fig['data'] if t.get('x')][0]
        self.assertIn('Casein beta chain', xs)
        self.assertNotIn('Beta-casein', xs)

    def test_rename_does_not_recolour_or_reorder(self):
        common = dict(selected_groups=['HighDose', 'Ctrl', 'LowDose'],
                      selected_proteins=['P02666', 'P02662'],
                      orientation='By Protein', plot_type='Stacked Bar Plots')
        base, _ = _plot(**common)
        renamed, _ = _plot(label_overrides={'sample_groups': {'HighDose': 'High'}},
                           **common)

        def colors(fig):
            return [(t.get('marker', {}) or {}).get('color') for t in fig['data']]

        # colours are keyed by the original group name → unchanged by a rename
        self.assertEqual(colors(base), colors(renamed))

        # the renamed group keeps its position in the trace order
        base_names = [t.get('name') for t in base['data']]
        new_names = [t.get('name') for t in renamed['data']]
        self.assertEqual(base_names.index('HighDose'), new_names.index('High'))
        self.assertNotIn('HighDose', new_names)


class TestSelectorOptionsCarryBothNames(unittest.TestCase):
    """The "Strip protein name" toggle lives next to the Proteins picker now and
    previews both forms client-side, so the options payload must ship both."""

    def test_protein_options_include_raw_and_stripped(self):
        df = pd.DataFrame({
            'Unique Peptide ID': ['p1', 'p2'],
            'Protein': ['P1', 'P2'],
            'protein_name': ['LACB_BOVIN Beta-lactoglobulin OS=Bos taurus GN=LGB',
                             'Beta-casein'],
            'Avg_A': [1.0, 2.0], 'Avg_B': [1.5, 2.5],
        })
        gdd, renamed, _w = dp.process_group_data(df)
        opts = dp.get_selector_options(renamed, gdd, dp.extract_protein_dict(renamed))
        by_id = {o['id']: o for o in opts['proteins']}
        self.assertEqual(by_id['P1']['name_raw'],
                         'LACB_BOVIN Beta-lactoglobulin OS=Bos taurus GN=LGB')
        self.assertEqual(by_id['P1']['name_stripped'], 'Beta-lactoglobulin')
        self.assertEqual(by_id['P2']['name_raw'], 'Beta-casein')
        self.assertEqual(by_id['P2']['name_stripped'], 'Beta-casein')


if __name__ == '__main__':
    unittest.main()
