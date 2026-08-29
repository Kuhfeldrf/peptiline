"""BLAST / homology search against the mbpdb_replica dataset.

The test harness can't stand up the replica DB (see the other DB-backed tests),
so these cover the query-*shape* fixes with a fake ORM: the exact-match pass
issues ONE ``filter(peptide__in=...)`` (not one query per peptide), no BLAST
FASTA / makeblastdb runs when nothing needs BLAST, and the per-search scratch
directory is always removed.
"""
import os
import shutil
import tempfile
import types
import unittest
from unittest import mock

from data_transformation.services import blast_search


class _FakeProtein:
    def __init__(self, pid, desc='desc', species='Bovine'):
        self.pid, self.desc, self.species = pid, desc, species


class _FakeRef:
    def __init__(self, doi):
        self.doi = doi
        self.additional_details = self.ptm = ''
        self.ic50 = self.inhibition_type = self.inhibited_microorganisms = None
        self.title = self.authors = self.abstract = ''


class _RelatedManager:
    def __init__(self, items):
        self._items = items

    def all(self):
        return list(self._items)


class _FakeFunction:
    def __init__(self, name, refs):
        self.function = name
        self.references = _RelatedManager(refs)


class _FakePeptide:
    def __init__(self, pk, seq, protein, functions):
        self.id = pk
        self.peptide = seq
        self.protein = protein
        self.intervals = '1-6'
        self.functions = _RelatedManager(functions)


class _FakeQS(list):
    """Minimal queryset: chainable no-ops + iterable of rows."""
    def filter(self, *a, **kw):
        return self

    def select_related(self, *a, **kw):
        return self

    def prefetch_related(self, *a, **kw):
        return self


def _fake_peptideinfo(rows, filter_spy=None):
    prot = _FakeProtein('P1')
    peps = [_FakePeptide(i + 1, s, prot,
                         [_FakeFunction('ACE-inhibitory', [_FakeRef('10.1/a')])])
            for i, s in enumerate(rows)]

    class _Objects:
        def filter(self, *a, **kw):
            if filter_spy is not None:
                filter_spy(kw)
            wanted = set(kw.get('peptide__in', []))
            return _FakeQS([p for p in peps if p.peptide in wanted])

    return types.SimpleNamespace(objects=_Objects())


class TestExactMatchPass(unittest.TestCase):
    def test_single_batched_filter_for_all_exact_peptides(self):
        calls = []
        fake = _fake_peptideinfo(['ABCDEF', 'GHIK', 'LMNOP'], filter_spy=calls.append)
        with mock.patch.object(blast_search, 'PeptideInfo', fake):
            df = blast_search.run_blast_search(['ABCDEF', 'GHIK', 'LMNOP', 'NOPE'],
                                               similarity_threshold=100)
        # ONE filter call, keyed by peptide__in — not one per sequence.
        self.assertEqual(len(calls), 1)
        self.assertIn('peptide__in', calls[0])
        self.assertEqual(set(df['search_peptide']), {'ABCDEF', 'GHIK', 'LMNOP'})

    def test_no_blast_db_built_for_exact_only_search(self):
        fake = _fake_peptideinfo(['ABCDEF'])
        work_dir = tempfile.mkdtemp(prefix='work_')
        try:
            with mock.patch.object(blast_search, 'PeptideInfo', fake), \
                 mock.patch.object(blast_search, 'make_blast_db') as mk:
                blast_search.run_blast_search(['ABCDEF'], similarity_threshold=100,
                                              work_dir=work_dir)
            mk.assert_not_called()
            self.assertEqual(os.listdir(work_dir), [])   # no scratch left behind
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)

    def test_empty_result_frame_keeps_columns(self):
        fake = _fake_peptideinfo([])
        with mock.patch.object(blast_search, 'PeptideInfo', fake):
            df = blast_search.run_blast_search(['ZZZ'], similarity_threshold=100)
        self.assertTrue(df.empty)
        self.assertIn('search_peptide', df.columns)
        self.assertIn('function', df.columns)

    def test_empty_replica_raises_actionable_error_on_blast_search(self):
        # threshold < 100 -> BLAST pass. An empty replica must raise a clear
        # message (surfaced to the UI), not a cryptic makeblastdb failure.
        fake = _fake_peptideinfo([])
        fake.objects.values_list = lambda *a, **k: types.SimpleNamespace(iterator=lambda: iter([]))
        work_dir = tempfile.mkdtemp(prefix='work_')
        try:
            with mock.patch.object(blast_search, 'PeptideInfo', fake), \
                 mock.patch.object(blast_search, 'make_blast_db') as mk:
                with self.assertRaises(RuntimeError) as ctx:
                    blast_search.run_blast_search(['LONGPEPTIDESEQ'], similarity_threshold=80,
                                                  work_dir=work_dir)
            self.assertIn('mbpdb_replica', str(ctx.exception))
            mk.assert_not_called()
            self.assertEqual(os.listdir(work_dir), [])   # scratch still cleaned
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)


class TestPeptideRows(unittest.TestCase):
    def _pep(self, functions):
        return _FakePeptide(7, 'ABCDEF', _FakeProtein('P9'), functions)

    def test_function_with_no_reference_emits_nothing(self):
        rows = blast_search._peptide_rows(self._pep([_FakeFunction('Antioxidant', [])]), 'ABCDEF')
        self.assertEqual(rows, [])

    def test_no_functions_emits_single_null_row(self):
        rows = blast_search._peptide_rows(self._pep([]), 'ABCDEF')
        self.assertEqual(len(rows), 1)
        self.assertIsNone(rows[0]['function'])
        self.assertEqual(rows[0]['protein_id'], 'P9')

    def test_blast_detail_merged_in(self):
        detail = ['95.00', '1', '6', '1', '6', '1e-3', '6', '0', '0']
        rows = blast_search._peptide_rows(
            self._pep([_FakeFunction('ACE-inhibitory', [_FakeRef('10.1/x')])]),
            'ABCDEF', blast_detail=detail)
        self.assertEqual(rows[0]['% Alignment'], '95.00')
        self.assertEqual(rows[0]['Gap opens'], '0')

    def test_as_int(self):
        self.assertEqual(blast_search._as_int('5'), 5)
        self.assertIsNone(blast_search._as_int('x'))
        self.assertIsNone(blast_search._as_int(None))


if __name__ == '__main__':
    unittest.main()
