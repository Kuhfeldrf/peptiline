"""
Tests for the "arrange & rename" panel wiring in the Heat Map dashboard:

  * ``selected_var_keys`` / ``selected_proteins`` order drives row order;
  * ``label_overrides`` (keyed by the ORIGINAL protein id / sample-group name)
    rewrite the display fields only — ``combo_key`` / ``protein_id`` / ``var_key``
    stay original so row identity, lookups and colour scales are untouched.

Run with:
    .venv/bin/python -m pytest tests/test_heatmap_labels.py -v
"""
import io
import unittest

import pandas as pd

from heatmap_viz.services import data_processor as hdp

SEQ = 'MKVLILACDEFGHIKLMNPQRSTVWY'


def _csv_bytes(df):
    buf = io.BytesIO()
    df.to_csv(buf, index=False)
    buf.seek(0)
    return buf


def _available(proteins, var_keys):
    df = pd.DataFrame({
        'Protein': ['P02666', 'P02666', 'P02662'],
        'protein_name': ['Beta-casein', 'Beta-casein', 'Alpha-S1-casein'],
        'protein_species': ['Bos taurus'] * 3,
        'start': [3, 10, 4],
        'end': [9, 16, 12],
        'function': ['bitter', None, None],
        'Unique Peptide ID': ['pep1', 'pep2', 'pep3'],
        'Avg_Bitter': [11.0, 21.0, 7.0],
        'Avg_NonBitter': [3.0, 5.0, 2.0],
    })
    merged, gdd, pdict, _col_order, err = hdp.load_merged_file(_csv_bytes(df), 'merged.csv')
    assert err is None, err
    pdict['P02666']['sequence'] = SEQ
    pdict['P02662']['sequence'] = SEQ
    avail, _msgs = hdp.build_available_data_variables(merged, pdict, gdd, proteins, var_keys)
    return avail


class TestRowOrder(unittest.TestCase):
    def test_var_key_order_drives_row_order(self):
        avail = _available(['P02666'], ['NonBitter', 'Bitter'])
        self.assertEqual([vd['var_label'] for vd in avail.values()],
                         ['NonBitter', 'Bitter'])


class TestApplyLabelOverrides(unittest.TestCase):
    def test_rename_display_fields_only(self):
        avail = _available(['P02666'], ['Bitter', 'NonBitter'])
        keys_before = list(avail.keys())

        hdp._apply_label_overrides(avail, {
            'sample_groups': {'Bitter': 'Bitter fraction'},
            'proteins': {'P02666': 'β-casein'},
        })

        # dict keys / identity fields untouched
        self.assertEqual(list(avail.keys()), keys_before)
        self.assertTrue(all(vd['protein_id'] == 'P02666' for vd in avail.values()))

        by_var = {k.split('_', 1)[1]: v for k, v in avail.items()}
        # single protein → label is the bare (renamed) sample name, no prefix
        self.assertEqual(by_var['Bitter']['var_label'], 'Bitter fraction')
        self.assertEqual(by_var['Bitter']['label'], 'Bitter fraction')
        self.assertEqual(by_var['NonBitter']['var_label'], 'NonBitter')  # not renamed
        self.assertEqual(by_var['Bitter']['protein_name'], 'β-casein')

    def test_multi_protein_label_uses_renamed_protein_prefix(self):
        avail = _available(['P02666', 'P02662'], ['Bitter'])
        hdp._apply_label_overrides(avail, {'proteins': {'P02662': 'αS1-casein'}})
        labels = {vd['label'] for vd in avail.values()}
        self.assertEqual(labels, {'Beta-casein Bitter', 'αS1-casein Bitter'})

    def test_none_overrides_is_noop(self):
        avail = _available(['P02666'], ['Bitter'])
        snapshot = {k: dict(v) for k, v in avail.items()}
        hdp._apply_label_overrides(avail, None)
        hdp._apply_label_overrides(avail, {})
        for k, v in avail.items():
            self.assertEqual(v['label'], snapshot[k]['label'])
            self.assertEqual(v['var_label'], snapshot[k]['var_label'])


if __name__ == '__main__':
    unittest.main()
