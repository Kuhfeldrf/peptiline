"""
Tests for the heatmap "Legend Position" option.

'right' keeps the historical single combined legend beside the plot; 'below'
splits each legend group into its own unit under the x-axis, the units laid
out side by side (wrapping to further rows when they would not fit).

Run with:
    cd /home/kuhfeldrf/mbpdb/include/peptide
    /home/kuhfeldrf/mbpdb/.venv/bin/python -m pytest tests/test_heatmap_legend_position.py -v
"""
import unittest

import numpy as np
import pandas as pd

# conftest.py bootstraps Django; import the service directly.
from heatmap_viz.services import heatmap_renderer as R

SEQ = 'MKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPF'
LEGEND_TITLES = ['Sample Type:', 'Peptide Counts:', 'Average Abundance:', '',
                 'Average Abundance:']
Y_TICKS = [1e4, 5e5, 1e6]


def _make_var(label, seed):
    rng = np.random.default_rng(seed)
    n = len(SEQ)
    counts = rng.integers(0, 8, n)
    ms = rng.random(n) * 1e6 + 1e4
    return {
        'protein_id': 'P02666', 'protein_name': 'Beta-casein', 'label': label,
        'var_label': label, 'protein_sequence': SEQ, 'AA_list': list(SEQ),
        'peptide_counts': pd.Series(counts), 'ms_data': pd.Series(ms),
        'ms_data_list': list(ms),
        'heatmap_df': pd.DataFrame({'count': counts, 'average': ms, 'AA': list(SEQ)}),
        'bioactive_peptide_abs_df': pd.DataFrame({'average': ms}),
        'bioactive_peptide_func_df': pd.DataFrame(),
    }


def _vars():
    return {'High': _make_var('High NE', 1), 'Low': _make_var('Low NE', 2)}


def _interactive(orientation, legend_position):
    import matplotlib.pyplot as plt
    R.max_count, R.num_colors = 8, 6
    return R.visualize_sequence_heatmap_interactive(
        _vars(), orientation, 'Beta-casein Sequence', 'Averaged Peptide Abundance',
        LEGEND_TITLES, plt.get_cmap('RdYlGn_r'), plt.get_cmap('Dark2'), 'Set3',
        [], [], 'yes', 'no', 'all-peptides', False, Y_TICKS,
        legend_position=legend_position,
    )


def _compact(legend_position):
    import matplotlib.pyplot as plt
    R.max_count, R.num_colors = 8, 6
    return R.visualize_sequence_heatmap_compact(
        _vars(), 'Beta-casein Sequence', 'Abundance', LEGEND_TITLES,
        plt.get_cmap('RdYlGn_r'), Y_TICKS, False,
        legend_position=legend_position,
    )


def _legends(fig):
    """All populated legend objects on the figure, keyed by layout name."""
    out = {}
    for key in ('legend', 'legend2', 'legend3', 'legend4'):
        leg = getattr(fig.layout, key, None)
        if leg is not None and leg.to_plotly_json():
            out[key] = leg
    return out


class TestNormalizeLegendPosition(unittest.TestCase):
    def test_below_aliases(self):
        for value in ('below', 'Below', ' bottom ', 'under'):
            self.assertEqual(R.normalize_legend_position(value), 'below')

    def test_default_and_unknown_fall_back_to_right(self):
        for value in (None, '', 'right', 'left', 'nonsense'):
            self.assertEqual(R.normalize_legend_position(value), 'right')


class TestRightIsUnchanged(unittest.TestCase):
    """'right' must keep the one combined, vertically-stacked legend."""

    def test_single_vertical_legend_with_group_titles(self):
        for orientation in ('Landscape', 'Portrait'):
            with self.subTest(orientation=orientation):
                fig = _interactive(orientation, 'right')
                legends = _legends(fig)
                self.assertEqual(list(legends), ['legend'])
                self.assertEqual(legends['legend'].orientation, 'v')
                # every visible entry carries its group header inside that legend
                titles = {t.legendgrouptitle.text for t in fig.data if t.showlegend}
                self.assertIn('Sample Type:', titles)
                self.assertIn('Peptide Counts:', titles)
                # nothing is routed to a secondary legend
                self.assertFalse([t for t in fig.data if t.legend not in (None, 'legend')])
                # width is still reserved on the right for it
                self.assertGreaterEqual(fig.layout.margin.r, 180)


class TestBelowSplitsIntoUnits(unittest.TestCase):
    def test_landscape_two_units_side_by_side(self):
        fig = _interactive('Landscape', 'below')
        legends = _legends(fig)
        self.assertEqual(sorted(legends), ['legend', 'legend2'])
        self.assertEqual([legends[k].title.text for k in ('legend', 'legend2')],
                         ['Sample Type:', 'Peptide Counts:'])
        for leg in legends.values():
            self.assertEqual(leg.orientation, 'h')     # items run across, not down
            self.assertLess(leg.y, 0)                  # under the x-axis
            self.assertEqual(leg.yanchor, 'top')
        # one row → same y, drawn left to right within the figure
        self.assertEqual(legends['legend'].y, legends['legend2'].y)
        self.assertLess(legends['legend'].x, legends['legend2'].x)
        for leg in legends.values():
            self.assertTrue(0 < leg.x < 1)
            self.assertEqual(leg.xanchor, 'center')

    def test_portrait_adds_the_abundance_unit(self):
        fig = _interactive('Portrait', 'below')
        titles = [leg.title.text for leg in _legends(fig).values()]
        self.assertEqual(sorted(titles),
                         ['Average Abundance:', 'Peptide Counts:', 'Sample Type:'])

    def test_compact_two_units(self):
        fig = _compact('below')
        titles = [leg.title.text for leg in _legends(fig).values()]
        self.assertEqual(sorted(titles), ['Average Abundance:', 'Peptide Counts:'])

    def test_entries_are_routed_to_their_unit(self):
        fig = _interactive('Landscape', 'below')
        by_legend = {}
        for trace in fig.data:
            if trace.showlegend:
                by_legend.setdefault(trace.legend, []).append(trace.name)
        self.assertEqual(sorted(by_legend['legend']), ['High NE', 'Low NE'])
        self.assertTrue(all('-' in n for n in by_legend['legend2']))
        # no legendgroup survives: inside a horizontal legend Plotly would lay
        # each group out as its own vertical column.
        self.assertFalse([t for t in fig.data if t.showlegend and t.legendgroup])

    def test_width_moves_from_the_right_margin_to_the_bottom(self):
        right = _interactive('Landscape', 'right')
        below = _interactive('Landscape', 'below')
        self.assertLess(below.layout.margin.r, right.layout.margin.r)
        self.assertGreater(below.layout.margin.b, right.layout.margin.b)
        # the plot rows themselves are untouched — the figure just grows taller
        self.assertGreater(below.layout.height, right.layout.height)


class TestUnitPacking(unittest.TestCase):
    """Row packing and the pixel→paper-fraction arithmetic behind the strip."""

    @staticmethod
    def _plan(n_units, n_entries=6):
        units = {
            f'Group {i}:': {'ref': 'legend' if i == 0 else f'legend{i + 1}',
                            'labels': [f'entry {j}' for j in range(n_entries)]}
            for i in range(n_units)
        }
        return R._below_legend_plan(units, 15)

    def test_units_wrap_to_a_second_row_when_too_wide(self):
        plan = self._plan(4)
        rows = R._pack_legend_rows(plan, R.LEGEND_REF_WIDTH_PX)
        self.assertGreater(len(rows), 1)
        self.assertEqual(sum(len(r) for r in rows), 4)
        for row in rows:
            span = (sum(u['width'] for u in row)
                    + R.LEGEND_UNIT_GAP_PX * (len(row) - 1))
            self.assertLessEqual(span, R.LEGEND_REF_WIDTH_PX)

    def test_stacked_puts_one_unit_per_row(self):
        plan = self._plan(3, n_entries=1)      # narrow units that WOULD share a row
        self.assertEqual(len(R._pack_legend_rows(plan, R.LEGEND_REF_WIDTH_PX)), 1)
        stacked = R._pack_legend_rows(plan, R.LEGEND_REF_WIDTH_PX, stacked=True)
        self.assertEqual([len(r) for r in stacked], [1, 1, 1])

    def test_stacked_never_crowds_however_narrow_the_figure(self):
        plan = self._plan(3)
        for width in (1600, 1000, 700, 420):
            positions, _h, _b, rows = R._below_legend_geometry(
                plan, True, width, top_px=30, base_height_px=500)
            self.assertEqual(rows, 3)
            # one unit per row → every unit centred, none beside another
            self.assertEqual({round(x, 6) for x, _y in positions.values()}, {0.5})
            self.assertEqual(len({y for _x, y in positions.values()}), 3)

    def test_rows_stack_downwards_below_the_axis(self):
        plan = self._plan(4)
        positions, height, bottom, rows = R._below_legend_geometry(
            plan, False, R.LEGEND_REF_WIDTH_PX, top_px=30, base_height_px=500)
        ys = sorted({y for _x, y in positions.values()}, reverse=True)
        self.assertEqual(len(ys), rows)
        self.assertTrue(all(y < 0 for y in ys))          # all under the x-axis
        self.assertEqual(height, 500 + R.LEGEND_ROW_PX * rows)
        self.assertEqual(bottom, R.LEGEND_BASE_BOTTOM_PX + R.LEGEND_ROW_PX * rows)

    def test_each_row_is_centred_on_the_figure(self):
        plan = self._plan(3)
        width = 1000
        positions, _h, _b, _rows = R._below_legend_geometry(
            plan, False, width, top_px=30, base_height_px=500)
        for row in R._pack_legend_rows(plan, width):
            # left edge of the first unit and right edge of the last, in the
            # same paper fractions Plotly positions them with
            left = positions[row[0]['ref']][0] - row[0]['width'] / (2 * width)
            right = positions[row[-1]['ref']][0] + row[-1]['width'] / (2 * width)
            self.assertAlmostEqual((left + right) / 2, 0.5, places=6)

    def test_units_never_overlap_at_the_width_they_were_packed_for(self):
        plan = self._plan(3)
        for width in (1600, 1200, 1000, 800, 600):
            positions, _h, _b, _rows = R._below_legend_geometry(
                plan, False, width, top_px=30, base_height_px=500)
            for row in R._pack_legend_rows(plan, width):
                edges = [(positions[u['ref']][0] * width - u['width'] / 2,
                          positions[u['ref']][0] * width + u['width'] / 2)
                         for u in row]
                for (_l1, r1), (l2, _r2) in zip(edges, edges[1:]):
                    self.assertGreaterEqual(l2 - r1, 0, f'overlap at width {width}')


class TestRepackForExportWidth(unittest.TestCase):
    """The static-export path re-lays the units for the width it renders at."""

    @staticmethod
    def _json(fig):
        import plotly.io as pio
        return pio.to_json(fig)

    def test_narrow_export_is_repacked_without_overlap(self):
        import json as _json
        fig = _interactive('Portrait', 'below')
        packed = R.repack_below_legends(self._json(fig), 700)
        layout = _json.loads(packed)['layout']
        plan = layout['meta']['below_legend']
        avail = 700 - layout['margin']['l'] - layout['margin']['r']

        by_y = {}
        for unit in plan['units']:
            leg = layout[unit['ref']]
            by_y.setdefault(round(leg['y'], 9), []).append((leg['x'], unit['width']))
        for row in by_y.values():
            row.sort()
            edges = [(x * avail - w / 2, x * avail + w / 2) for x, w in row]
            for (_l1, r1), (l2, _r2) in zip(edges, edges[1:]):
                self.assertGreaterEqual(l2 - r1, 0)
        # height and bottom margin follow the row count it settled on
        self.assertEqual(layout['height'],
                         plan['base_height_px'] + R.LEGEND_ROW_PX * len(by_y))
        self.assertEqual(layout['margin']['b'],
                         plan['base_bottom_px'] + R.LEGEND_ROW_PX * len(by_y))

    def test_right_legend_figures_pass_through_untouched(self):
        original = self._json(_interactive('Landscape', 'right'))
        self.assertEqual(R.repack_below_legends(original, 700), original)

    def test_missing_width_or_bad_json_is_a_no_op(self):
        original = self._json(_interactive('Landscape', 'below'))
        self.assertEqual(R.repack_below_legends(original, None), original)
        self.assertEqual(R.repack_below_legends('not json', 700), 'not json')


class TestStackedFigures(unittest.TestCase):
    def test_stacked_gives_every_unit_its_own_row(self):
        fig = _interactive('Portrait', 'below-stacked')
        legends = _legends(fig)
        self.assertEqual(len(legends), 3)
        self.assertEqual({round(leg.x, 6) for leg in legends.values()}, {0.5})
        self.assertEqual(len({leg.y for leg in legends.values()}), 3)

    def test_stacked_is_taller_than_side_by_side(self):
        self.assertGreater(_interactive('Portrait', 'below-stacked').layout.height,
                           _interactive('Portrait', 'below').layout.height)

    def test_compact_and_landscape_accept_stacked(self):
        for fig in (_compact('below-stacked'), _interactive('Landscape', 'below-stacked')):
            legends = _legends(fig)
            self.assertEqual({round(leg.x, 6) for leg in legends.values()}, {0.5})


if __name__ == '__main__':
    unittest.main()
