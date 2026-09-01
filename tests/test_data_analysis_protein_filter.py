"""
Regression: a "Selected Protein(s)" plot must show ONLY the proteins the user
picked. filter_dataframe() keeps a peptide row when any of its semicolon-listed
proteins is selected, so a peptide shared between a selected and an unselected
protein used to drag the unselected partner onto the figure (casein isoforms
share many peptides). process_protein_data() now restricts protein_df to the
selection unless plot_minor / "All Proteins" is in play.

Run with:
    .venv/bin/python -m pytest tests/test_data_analysis_protein_filter.py -v
"""
import json
import unittest

import numpy as np
import pandas as pd

from data_analysis.services import plotter

GROUPS = ['Ctrl', 'Treat']


def _dataset():
    """P1 & P2 share pep0; pep1 is P1-only, pep2 is P2-only, pep3 is P3-only."""
    rng = np.random.default_rng(1)
    reps = {g: [f'{g}_{i}' for i in range(1, 4)] for g in GROUPS}
    spec = [('pep0', 'P1;P2'), ('pep1', 'P1'), ('pep2', 'P2'), ('pep3', 'P3')]
    rows = []
    for pep, prot in spec:
        row = {'Unique Peptide ID': pep, 'Protein': prot}
        for g in GROUPS:
            for rc in reps[g]:
                row[rc] = float(100 + rng.normal(0, 5))
            row[f'Avg_{g}'] = float(np.mean([row[rc] for rc in reps[g]]))
        rows.append(row)
    protein_dict = {'P1': {'name': 'Beta-casein'},
                    'P2': {'name': 'Alpha-S1-casein'},
                    'P3': {'name': 'Kappa-casein'}}
    return pd.DataFrame(rows), reps, protein_dict


def _bar_categories(fig):
    """Every distinct category label drawn on the x-axis / bar hover."""
    names = set()
    for t in fig['data']:
        for h in (t.get('hovertext') or []):
            if isinstance(h, str) and h.startswith('Protein: '):
                names.add(h.split('Protein: ', 1)[1].split('<br>', 1)[0])
    return names


def _plot(**extra):
    merged, gdd, pdict = _dataset()
    params = dict(
        selected_groups=GROUPS,
        selected_functions=['No Function Filter'],
        plot_type='Grouped Bar Plots', orientation='By Protein',
        abs_or_count='Abundance', metric_type='Absolute',
    )
    params.update(extra)
    fig_json, warns = plotter.generate_plot(merged, gdd, pdict, params)
    assert fig_json, warns
    return json.loads(fig_json)


class TestSelectedProteinsOnly(unittest.TestCase):
    def test_shared_peptide_does_not_leak_unselected_protein(self):
        fig = _plot(selected_proteins=['P1'])
        cats = _bar_categories(fig)
        self.assertIn('Beta-casein', cats)
        self.assertNotIn('Alpha-S1-casein', cats)   # shared pep0 partner
        self.assertNotIn('Kappa-casein', cats)

    def test_two_selected_no_third(self):
        fig = _plot(selected_proteins=['P1', 'P2'])
        cats = _bar_categories(fig)
        self.assertEqual(cats, {'Beta-casein', 'Alpha-S1-casein'})

    def test_plot_minor_still_pools_the_rest(self):
        fig = _plot(selected_proteins=['P1'], plot_minor=True)
        cats = _bar_categories(fig)
        self.assertIn('Beta-casein', cats)
        self.assertIn('Other Proteins', cats)
        self.assertNotIn('Alpha-S1-casein', cats)


if __name__ == '__main__':
    unittest.main()
