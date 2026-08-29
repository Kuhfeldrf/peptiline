"""
Standalone tests for the data_transformation service modules.

These tests run without a live Django instance by mocking settings.
They verify that the Django services produce output matching the logic
from the original Jupyter notebook (data_transformation.ipynb).

Run with:
    cd /home/kuhfeldrf/mbpdb/include/peptide
    /home/kuhfeldrf/mbpdb/.venv/bin/python -m pytest tests/test_data_transformation.py -v
"""
"""
Standalone tests for the data_transformation service modules.

These tests run without a live Django instance.  Django is bootstrapped
via tests/conftest.py (sets up settings.configure() and blocks Celery).

Run with:
    cd /home/kuhfeldrf/mbpdb/include/peptide
    /home/kuhfeldrf/mbpdb/.venv/bin/python -m pytest tests/test_data_transformation.py -v
"""
import io
import json
import sys
import os
import unittest

import pandas as pd
import numpy as np

# conftest.py already bootstrapped Django — just import the services.
from data_transformation.services import data_loader, group_processing, data_combiner, export_manager, protein_handler


# ---------------------------------------------------------------------------
# Helpers – synthetic test data
# ---------------------------------------------------------------------------

def make_peptidomic_df():
    """Synthetic peptidomic DataFrame mimicking Proteome Discoverer output."""
    return pd.DataFrame({
        'Annotated Sequence': [
            '[A].RELEELNVPGEI.[V]',
            '[S].SEESITRINK.[K]',
            '[A].RELEELNVPGEIVES.[L]',
            '[N].SLPQNIPPLTQ.[G]',
        ],
        'Modifications': [
            None,
            '1xPhospho [S1]',
            '1xPhospho [S15]',
            None,
        ],
        'Master Protein Accessions': [
            'P02666',
            'P02666',
            'P02666',
            'Q9NR61',
        ],
        'Positions in Proteins': [
            'P02666 [16-27]',
            'P02666 [34-43]',
            'P02666 [16-30]',
            'Q9NR61 [84-94]',
        ],
        'Sequence': [
            'RELEELNVPGEI',
            'SEESITRINK',
            'RELEELNVPGEIVES',
            'SLPQNIPPLTQ',
        ],
        'Sample_A1': [2257142.857, 1816666.667, 3083333.333, 12645000.0],
        'Sample_A2': [1900000.0,  1500000.0,   2800000.0,   10000000.0],
        'Sample_B1': [800000.0,   600000.0,    1200000.0,   5000000.0],
        'Sample_B2': [700000.0,   500000.0,    1000000.0,   4500000.0],
    })


def make_mbpdb_df():
    """Synthetic MBPDB search results DataFrame."""
    return pd.DataFrame({
        'search_peptide': ['RELEELNVPGEI', 'SEESITRINK', 'RELEELNVPGEI'],
        'peptide':        ['RELEELNVPGEI', 'SEESITRINK', 'RELEELNVPGEI'],
        'protein_id':     ['P02666',       'P02666',     'P02666'],
        'protein_description': ['Beta-casein', 'Beta-casein', 'Beta-casein'],
        'species':        ['Bovine',       'Bovine',     'Bovine'],
        'intervals':      ['[16-27]',      '[34-43]',    '[16-27]'],
        'function':       ['Opioid',       'ACE inhibitor', 'Opioid'],
        '% Alignment':    [100.0,          100.0,        100.0],
    })


def make_group_data():
    return {
        'GroupA': {'grouping_variable': 'GroupA', 'abundance_columns': ['Sample_A1', 'Sample_A2']},
        'GroupB': {'grouping_variable': 'GroupB', 'abundance_columns': ['Sample_B1', 'Sample_B2']},
    }


def make_protein_dict():
    return {
        'P02666': {'name': 'CASB_BOVIN Beta-casein', 'species': 'Bovine'},
        'Q9NR61': {'name': 'SOME_HUMAN Protein', 'species': 'Human'},
    }


# ---------------------------------------------------------------------------
# data_loader tests
# ---------------------------------------------------------------------------

class TestFindSpecies(unittest.TestCase):

    def test_bovine_detected(self):
        result = data_loader.find_species(">sp|P02666|CASB_BOVIN Beta-casein OS=Bos taurus")
        self.assertEqual(result, "Bovine")

    def test_human_detected(self):
        result = data_loader.find_species(">sp|P04406|G3P_HUMAN OS=Homo sapiens")
        self.assertEqual(result, "Human")

    def test_unknown_returned(self):
        result = data_loader.find_species(">sp|P12345|SOME_YEAST OS=Saccharomyces cerevisiae")
        self.assertEqual(result, "unknown")


class TestParseFastaHeaders(unittest.TestCase):

    FASTA = """>sp|P02666|CASB_BOVIN Beta-casein OS=Bos taurus OX=9913
RELEELNVPGEI
>sp|P04406|G3P_HUMAN Glyceraldehyde-3-phosphate dehydrogenase OS=Homo sapiens OX=9606
MVKVGVNG
"""

    def test_parse_uploaded_fasta(self):
        result = data_loader.parse_uploaded_fasta(self.FASTA.encode())
        self.assertIn('P02666', result)
        self.assertIn('P04406', result)
        self.assertEqual(result['P02666']['species'], 'Bovine')
        self.assertEqual(result['P04406']['species'], 'Human')
        self.assertIn('RELEELNVPGEI', result['P02666']['sequence'])

    def test_empty_fasta_returns_empty(self):
        result = data_loader.parse_uploaded_fasta(b'')
        self.assertEqual(result, {})


class TestValidateFastaFormat(unittest.TestCase):

    def test_valid_fasta(self):
        fasta = ">sp|P02666|CASB_BOVIN\nRELEELNVPGEI\n"
        valid, msg = data_loader.validate_fasta_format(fasta.encode())
        self.assertTrue(valid)
        self.assertEqual(msg, "")

    def test_no_header(self):
        valid, msg = data_loader.validate_fasta_format(b"RELEELNVPGEI\n")
        self.assertFalse(valid)

    def test_empty_file(self):
        valid, msg = data_loader.validate_fasta_format(b"")
        self.assertFalse(valid)


class TestExtractProteinId(unittest.TestCase):

    def test_uniprot_pipe_format(self):
        result = data_loader.extract_protein_id("sp|P02666|CASB_BOVIN")
        self.assertEqual(result, "P02666")

    def test_peaks_two_part_pipe(self):
        """PEAKS uses ACCESSION|NAME (2 fields); accession is the first field."""
        result = data_loader.extract_protein_id("A0A087WWV8|A0A087WWV8_HUMAN")
        self.assertEqual(result, "A0A087WWV8")

    def test_peaks_semicolon_pipe_multi(self):
        """PEAKS lists multiple proteins with ';' between ACCESSION|NAME pairs."""
        result = data_loader.extract_protein_id(
            "A0A087WWV8|A0A087WWV8_HUMAN;P01834|IGKC_HUMAN"
        )
        self.assertEqual(result, ["A0A087WWV8", "P01834"])

    def test_plain_accession(self):
        result = data_loader.extract_protein_id("P02666")
        self.assertEqual(result, "P02666")

    def test_semicolon_separated(self):
        result = data_loader.extract_protein_id("P02666; Q9NR61")
        self.assertIsInstance(result, list)
        self.assertIn("P02666", result)
        self.assertIn("Q9NR61", result)

    def test_nan_returns_empty(self):
        result = data_loader.extract_protein_id(float('nan'))
        self.assertEqual(result, [])

    def test_empty_string_returns_empty(self):
        result = data_loader.extract_protein_id('')
        self.assertEqual(result, [])


class TestLoadAndValidateFile(unittest.TestCase):

    def _make_csv_bytes(self, df):
        buf = io.BytesIO()
        df.to_csv(buf, index=False)
        return buf.getvalue()

    def test_valid_peptidomic_csv(self):
        df = make_peptidomic_df()
        content = self._make_csv_bytes(df)
        result_df, status, err, warn = data_loader.load_and_validate_file(
            content, 'test.csv', 'Peptidomic'
        )
        self.assertEqual(status, 'yes')
        self.assertIsNotNone(result_df)

    def test_valid_mbpdb_csv(self):
        df = make_mbpdb_df()
        content = self._make_csv_bytes(df)
        result_df, status, err, warn = data_loader.load_and_validate_file(
            content, 'test.csv', 'MBPDB'
        )
        self.assertEqual(status, 'yes')
        self.assertIsNotNone(result_df)

    def test_mbpdb_missing_required_column(self):
        df = pd.DataFrame({'search_peptide': ['ABC'], 'peptide': ['ABC']})
        content = self._make_csv_bytes(df)
        _, status, err, _ = data_loader.load_and_validate_file(
            content, 'test.csv', 'MBPDB'
        )
        self.assertEqual(status, 'no')
        self.assertIn('function', err)

    def test_unsupported_extension(self):
        _, status, err, _ = data_loader.load_and_validate_file(
            b'data', 'test.xyz', 'Peptidomic'
        )
        self.assertEqual(status, 'no')

    def test_column_mapping_applied(self):
        """Master Protein Accessions should be remapped to Protein."""
        df = make_peptidomic_df()
        content = self._make_csv_bytes(df)
        result_df, status, _, _ = data_loader.load_and_validate_file(
            content, 'test.csv', 'Peptidomic'
        )
        self.assertEqual(status, 'yes')
        self.assertIn('Protein', result_df.columns)

    def test_peaks_accession_column_mapped_and_split(self):
        """PEAKS 'Accession' column (';'-joined ACCESSION|NAME pairs) maps to
        Protein and splits into clean accessions; inline mods stripped."""
        df = pd.DataFrame({
            'Peptide': ['ADYEKHKVYAC(+57.02)EVTHQG', 'AENTLQSFR'],
            'Accession': [
                'A0A087WWV8|A0A087WWV8_HUMAN;P01834|IGKC_HUMAN',
                'P02666|CASB_BOVIN',
            ],
            'Area': [100, 200],
        })
        content = self._make_csv_bytes(df)
        result_df, status, _, _ = data_loader.load_and_validate_file(
            content, 'peaks.csv', 'Peptidomic'
        )
        self.assertEqual(status, 'yes')
        self.assertIn('Protein', result_df.columns)
        # multi-protein row split into clean accessions
        self.assertEqual(result_df['Protein'].iloc[0], ['A0A087WWV8', 'P01834'])
        self.assertEqual(result_df['Protein'].iloc[1], 'P02666')
        # inline (+57.02) modification stripped from the sequence
        self.assertEqual(result_df['Sequence'].iloc[0], 'ADYEKHKVYACEVTHQG')

    def test_literal_protein_column_normalized(self):
        """A column literally named 'Protein' with full 'sp|ACC|NAME'
        accessions (Skyline / native PEPEX) is reduced to bare IDs."""
        df = pd.DataFrame({
            'Peptide': ['YLLPEVTVLDYGKK'],
            'Protein': ['sp|O15194|CTDSL_HUMAN'],
            'Area': [100],
        })
        content = self._make_csv_bytes(df)
        result_df, status, _, _ = data_loader.load_and_validate_file(
            content, 'skyline.csv', 'Peptidomic'
        )
        self.assertEqual(status, 'yes')
        self.assertEqual(result_df['Protein'].iloc[0], 'O15194')


class TestDetectLongFormat(unittest.TestCase):
    """PEAKS/Spectronaut-style long (one row per peptide per run) exports
    must be recognised; already-wide exports (MaxQuant/Skyline/PEPEX-style)
    must not be, even when the same base sequence repeats across rows."""

    def test_peaks_style_long_data_detected(self):
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK', 'PEPTIDEK', 'SAMPLER', 'SAMPLER'],
            'Accession': ['A0A087WWV8|A0A087WWV8_HUMAN'] * 4,
            'Area': [100, 200, 300, 400],
            'Source File': ['run1.d', 'run2.d', 'run1.d', 'run2.d'],
        })
        needs_pivot, sample_col, value_col = data_loader.detect_long_format(df)
        self.assertTrue(needs_pivot)
        self.assertEqual(sample_col, 'Source File')
        self.assertEqual(value_col, 'Area')

    def test_spectronaut_style_long_data_detected(self):
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK', 'PEPTIDEK', 'SAMPLER', 'SAMPLER'],
            'EG.ModifiedPeptide': ['_PEPTIDEK_'] * 2 + ['_SAMPLER_'] * 2,
            'PG.ProteinGroups': ['P02666'] * 4,
            'EG.TotalQuantity (Settings)': [10.0, 20.0, 30.0, 40.0],
            'R.FileName': ['run1', 'run2', 'run1', 'run2'],
        })
        needs_pivot, sample_col, value_col = data_loader.detect_long_format(df)
        self.assertTrue(needs_pivot)
        self.assertEqual(sample_col, 'R.FileName')
        self.assertEqual(value_col, 'EG.TotalQuantity (Settings)')

    def test_wide_maxquant_style_data_not_detected(self):
        df = make_peptidomic_df().rename(columns={'Master Protein Accessions': 'Proteins'})
        needs_pivot, sample_col, value_col = data_loader.detect_long_format(df)
        self.assertFalse(needs_pivot)
        self.assertIsNone(sample_col)
        self.assertIsNone(value_col)

    def test_no_sample_column_not_detected(self):
        """Area/Intensity columns alone (no run identifier) must not trigger a pivot."""
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK', 'PEPTIDEK'],
            'Area': [100, 200],
        })
        needs_pivot, _, _ = data_loader.detect_long_format(df)
        self.assertFalse(needs_pivot)

    def test_skyline_style_repeated_base_sequence_not_detected(self):
        """Different charge states/mods of the same base 'Peptide' share a
        Source File-like column only coincidentally — since the Unique
        Peptide ID (from 'Peptide Modified Sequence') is already one row per
        peptide, this must not be treated as long format."""
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK', 'PEPTIDEK'],
            'Peptide Modified Sequence': ['PEPTIDEK', 'PEPTIDEK[+16]'],
            'Source File': ['run1.d', 'run1.d'],
            'Area': [100, 200],
        })
        needs_pivot, _, _ = data_loader.detect_long_format(df)
        self.assertFalse(needs_pivot)


class TestPivotLongToWide(unittest.TestCase):

    def test_peaks_style_pivot_sums_psms_within_sample(self):
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK', 'PEPTIDEK', 'PEPTIDEK', 'SAMPLER'],
            'Protein': ['P02666', 'P02666', 'P02666', 'Q9NR61'],
            'Area': [100, 50, 300, 400],
            'Source File': ['run1.d', 'run1.d', 'run2.d', 'run1.d'],
        })
        result = data_loader.pivot_long_to_wide(df, 'Source File', 'Area')
        self.assertEqual(len(result), 2)
        row = result[result['Sequence'] == 'PEPTIDEK'].iloc[0]
        # Two PSMs of PEPTIDEK in run1.d (100 + 50) must be summed
        self.assertEqual(row['run1.d'], 150)
        self.assertEqual(row['run2.d'], 300)
        sampler_row = result[result['Sequence'] == 'SAMPLER'].iloc[0]
        self.assertTrue(pd.isna(sampler_row['run2.d']))

    def test_spectronaut_style_pivot_does_not_double_count_elution_group(self):
        """EG.TotalQuantity (Settings) is already reported once per
        (modified peptide, run) — residual duplicate rows (e.g. multiple
        charge states) must be collapsed with max, not summed."""
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK', 'PEPTIDEK'],
            'EG.ModifiedPeptide': ['_PEPTIDEK_', '_PEPTIDEK_'],
            'PG.ProteinGroups': ['P02666', 'P02666'],
            'EG.TotalQuantity (Settings)': [500.0, 500.0],
            'R.FileName': ['run1', 'run1'],
        })
        result = data_loader.pivot_long_to_wide(df, 'R.FileName', 'EG.TotalQuantity (Settings)')
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]['run1'], 500.0)

    def test_pivot_drops_non_schema_psm_columns(self):
        """Per-PSM descriptors with no Table 1 mapping (e.g. PEAKS '-10LgP',
        'Scan') must be dropped, not collapsed to an arbitrary first value
        that would then look like a spurious abundance column."""
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK', 'PEPTIDEK'],
            'Protein': ['P02666', 'P02666'],
            'Area': [100, 200],
            'Source File': ['run1.d', 'run2.d'],
            '-10LgP': [30.5, 31.2],
            'Scan': [111, 222],
        })
        result = data_loader.pivot_long_to_wide(df, 'Source File', 'Area')
        self.assertNotIn('-10LgP', result.columns)
        self.assertNotIn('Scan', result.columns)


class TestLoadAndValidateFileLongFormatPivot(unittest.TestCase):
    """End-to-end: uploading a long-format export auto-pivots before the
    rest of the pipeline sees it, and the resulting sample columns are
    correctly picked up as abundance columns."""

    def _make_csv_bytes(self, df):
        buf = io.BytesIO()
        df.to_csv(buf, index=False)
        return buf.getvalue()

    def test_peaks_long_export_pivots_and_is_flagged(self):
        df = pd.DataFrame({
            'Peptide': ['PEPTIDEK', 'PEPTIDEK', 'SAMPLER', 'SAMPLER'],
            'Accession': ['A0A087WWV8|A0A087WWV8_HUMAN'] * 4,
            'Area': [100, 200, 300, 400],
            'Source File': ['run1.d', 'run2.d', 'run1.d', 'run2.d'],
        })
        content = self._make_csv_bytes(df)
        result_df, status, err, warn = data_loader.load_and_validate_file(
            content, 'peaks.csv', 'Peptidomic'
        )
        self.assertEqual(status, 'yes', msg=err)
        self.assertIn('long-format', warn)
        self.assertEqual(len(result_df), 2)  # PEPTIDEK, SAMPLER
        filtered = data_loader.get_filtered_columns(result_df)
        self.assertEqual(set(filtered), {'run1.d', 'run2.d'})

    def test_spectronaut_long_export_pivots_and_is_flagged(self):
        df = pd.DataFrame({
            'PEP.StrippedSequence': ['PEPTIDEK', 'PEPTIDEK', 'SAMPLER', 'SAMPLER'],
            'EG.ModifiedPeptide': ['_PEPTIDEK_'] * 2 + ['_SAMPLER_'] * 2,
            'PG.ProteinGroups': ['P02666'] * 4,
            'EG.TotalQuantity (Settings)': [10.0, 20.0, 30.0, 40.0],
            'R.FileName': ['run1', 'run2', 'run1', 'run2'],
        })
        content = self._make_csv_bytes(df)
        result_df, status, err, warn = data_loader.load_and_validate_file(
            content, 'spectronaut.csv', 'Peptidomic'
        )
        self.assertEqual(status, 'yes', msg=err)
        self.assertIn('long-format', warn)
        self.assertEqual(len(result_df), 2)
        filtered = data_loader.get_filtered_columns(result_df)
        self.assertEqual(set(filtered), {'run1', 'run2'})

    def test_wide_export_not_pivoted(self):
        df = make_peptidomic_df().rename(columns={'Master Protein Accessions': 'Proteins'})
        content = self._make_csv_bytes(df)
        result_df, status, err, warn = data_loader.load_and_validate_file(
            content, 'wide.csv', 'Peptidomic'
        )
        self.assertEqual(status, 'yes', msg=err)
        self.assertNotIn('long-format', warn)
        self.assertEqual(len(result_df), len(df))


class TestGetFilteredColumnsAccessionAndId(unittest.TestCase):
    """Regression coverage for the exclusion-list fixes needed to support
    PEAKS long-format uploads: 'Accession' is a Table 1 protein-ID source
    column (always metadata), while the literal 'id' column (MaxQuant) must
    only be excluded on an exact match — not as a substring, which used to
    false-positive on sample names like a PEAKS Source File containing
    '...IDB001031...'."""

    def test_accession_excluded(self):
        df = pd.DataFrame({'Sequence': ['A'], 'Accession': ['P1|NAME'], 'Sample1': [100]})
        filtered = data_loader.get_filtered_columns(df)
        self.assertNotIn('Accession', filtered)
        self.assertIn('Sample1', filtered)

    def test_literal_id_column_excluded(self):
        df = pd.DataFrame({'Sequence': ['A'], 'id': [0], 'Sample1': [100]})
        filtered = data_loader.get_filtered_columns(df)
        self.assertNotIn('id', filtered)

    def test_sample_name_containing_id_substring_not_excluded(self):
        df = pd.DataFrame({
            'Sequence': ['A'],
            '20230208_TTS1_ptran0117_OBJ42352_IDB001031_9D1_Boco_IP1_L243_MAPPs_1_1_1419.d': [100],
        })
        filtered = data_loader.get_filtered_columns(df)
        self.assertEqual(len(filtered), 1)


class TestStripInlineModifications(unittest.TestCase):

    def test_peaks_parenthetical_mod(self):
        self.assertEqual(
            data_loader.strip_inline_modifications('ADYEKHKVYAC(+57.02)EVTHQG'),
            'ADYEKHKVYACEVTHQG',
        )

    def test_bracket_annotation(self):
        self.assertEqual(
            data_loader.strip_inline_modifications('NILREKQTDEIK[DN]'),
            'NILREKQTDEIK',
        )

    def test_enzymatic_flanks(self):
        self.assertEqual(
            data_loader.strip_inline_modifications('K.PEPTIDE.R'), 'PEPTIDE'
        )

    def test_clean_sequence_unchanged(self):
        self.assertEqual(
            data_loader.strip_inline_modifications('PEPTIDE'), 'PEPTIDE'
        )

    def test_non_string_passthrough(self):
        self.assertTrue(pd.isna(
            data_loader.strip_inline_modifications(float('nan'))
        ))


class TestExtractSequences(unittest.TestCase):

    def test_from_sequence_column(self):
        df = pd.DataFrame({'Sequence': ['ABCDE', 'FGHIJ', 'ABCDE']})
        seqs = data_loader.extract_sequences(df)
        self.assertIn('ABCDE', seqs)
        self.assertIn('FGHIJ', seqs)
        self.assertEqual(len(seqs), 2)  # unique only

    def test_from_annotated_sequence(self):
        df = pd.DataFrame({
            'Annotated Sequence': ['[A].HELLO.[V]', '[B].WORLD.[K]']
        })
        seqs = data_loader.extract_sequences(df)
        self.assertIn('HELLO', seqs)
        self.assertIn('WORLD', seqs)


class TestGetFilteredColumns(unittest.TestCase):

    def test_excludes_metadata_columns(self):
        df = make_peptidomic_df()
        filtered = data_loader.get_filtered_columns(df)
        self.assertNotIn('Sequence', filtered)
        self.assertNotIn('Master Protein Accessions', filtered)
        self.assertNotIn('Positions in Proteins', filtered)
        self.assertNotIn('Modifications', filtered)
        self.assertNotIn('Annotated Sequence', filtered)

    def test_retains_abundance_columns(self):
        df = make_peptidomic_df()
        filtered = data_loader.get_filtered_columns(df)
        self.assertIn('Sample_A1', filtered)
        self.assertIn('Sample_A2', filtered)
        self.assertIn('Sample_B1', filtered)
        self.assertIn('Sample_B2', filtered)


class TestFindMissingProteins(unittest.TestCase):

    def test_finds_missing(self):
        df = make_peptidomic_df()
        # Rename column to match validated format
        df = df.rename(columns={'Master Protein Accessions': 'Protein'})
        protein_dict = {'P02666': {'name': 'Beta-casein', 'species': 'Bovine'}}
        missing = data_loader.find_missing_proteins(df, protein_dict)
        self.assertIn('Q9NR61', missing)
        self.assertNotIn('P02666', missing)

    def test_empty_df_returns_empty(self):
        missing = data_loader.find_missing_proteins(pd.DataFrame(), {})
        self.assertEqual(missing, set())


# ---------------------------------------------------------------------------
# group_processing tests
# ---------------------------------------------------------------------------

class TestParseGroupJson(unittest.TestCase):

    def test_valid_json(self):
        group_json = json.dumps({
            'GroupA': ['Sample_A1', 'Sample_A2'],
            'GroupB': ['Sample_B1', 'Sample_B2'],
        }).encode()
        available = ['Sample_A1', 'Sample_A2', 'Sample_B1', 'Sample_B2']
        group_data, err = group_processing.parse_group_json(group_json, available)
        self.assertIsNone(err)
        self.assertEqual(len(group_data), 2)
        self.assertEqual(group_data['GroupA']['grouping_variable'], 'GroupA')
        self.assertEqual(group_data['GroupA']['abundance_columns'], ['Sample_A1', 'Sample_A2'])

    def test_missing_column_returns_error(self):
        group_json = json.dumps({'GroupA': ['Missing_Col']}).encode()
        _, err = group_processing.parse_group_json(group_json, ['Sample_A1'])
        self.assertIsNotNone(err)
        self.assertIn('Missing_Col', err)

    def test_invalid_json_returns_error(self):
        _, err = group_processing.parse_group_json(b'not valid json', [])
        self.assertIsNotNone(err)


class TestBuildGroupData(unittest.TestCase):

    def test_adds_group(self):
        gd, err = group_processing.build_group_data(
            ['Sample_A1', 'Sample_A2'], 'GroupA'
        )
        self.assertIsNone(err)
        self.assertEqual(gd['GroupA']['grouping_variable'], 'GroupA')
        self.assertEqual(gd['GroupA']['abundance_columns'], ['Sample_A1', 'Sample_A2'])

    def test_appends_to_existing(self):
        gd = {'GroupA': {'grouping_variable': 'GroupA', 'abundance_columns': ['Sample_A1']}}
        gd, err = group_processing.build_group_data(['Sample_B1'], 'GroupB', gd)
        self.assertIsNone(err)
        self.assertIn('GroupB', gd)
        self.assertEqual(gd['GroupB']['grouping_variable'], 'GroupB')

    def test_empty_name_returns_error(self):
        _, err = group_processing.build_group_data(['Sample_A1'], '')
        self.assertIsNotNone(err)

    def test_empty_columns_returns_error(self):
        _, err = group_processing.build_group_data([], 'GroupA')
        self.assertIsNotNone(err)


class TestBuildNoGroupData(unittest.TestCase):

    def test_creates_individual_groups(self):
        cols = ['Sample_A1', 'Sample_A2', 'Sample_B1']
        gd, err = group_processing.build_no_group_data(cols)
        self.assertIsNone(err)
        self.assertEqual(len(gd), 3)
        # Each column is its own group
        all_vars = {v['grouping_variable'] for v in gd.values()}
        self.assertEqual(all_vars, set(cols))

    def test_empty_columns_returns_error(self):
        _, err = group_processing.build_no_group_data([])
        self.assertIsNotNone(err)


class TestExportGroupDataJson(unittest.TestCase):

    def test_round_trip(self):
        gd = make_group_data()
        json_str = group_processing.export_group_data_json(gd)
        parsed = json.loads(json_str)
        self.assertIn('GroupA', parsed)
        self.assertEqual(parsed['GroupA'], ['Sample_A1', 'Sample_A2'])
        self.assertIn('GroupB', parsed)


# ---------------------------------------------------------------------------
# data_combiner tests
# ---------------------------------------------------------------------------

class TestAddProteinInfo(unittest.TestCase):

    def test_inserts_species_and_name(self):
        df = pd.DataFrame({
            'Protein': ['P02666', 'Q9NR61'],
            'Sequence': ['RELEELNVPGEI', 'SLPQNIPPLTQ'],
        })
        protein_dict = make_protein_dict()
        result = data_combiner.add_protein_info(df, protein_dict)
        self.assertIn('protein_species', result.columns)
        self.assertIn('protein_name', result.columns)
        self.assertEqual(result.loc[0, 'protein_species'], 'Bovine')
        self.assertEqual(result.loc[1, 'protein_species'], 'Human')

    def test_protein_columns_inserted_after_protein(self):
        df = pd.DataFrame({'Protein': ['P02666'], 'Sequence': ['PEPTIDE']})
        result = data_combiner.add_protein_info(df, make_protein_dict())
        cols = list(result.columns)
        protein_idx = cols.index('Protein')
        species_idx = cols.index('protein_species')
        name_idx = cols.index('protein_name')
        self.assertEqual(species_idx, protein_idx + 1)
        self.assertEqual(name_idx, protein_idx + 2)

    def test_unknown_protein_stays_unknown(self):
        df = pd.DataFrame({'Protein': ['UNKNOWN_ID'], 'Sequence': ['PEPTIDE']})
        result = data_combiner.add_protein_info(df, {})
        self.assertEqual(result.loc[0, 'protein_species'], 'Unknown')


class TestExtractBioactivePeptides(unittest.TestCase):

    def test_groups_by_search_peptide(self):
        mbpdb = make_mbpdb_df()
        cleaned, grouped = data_combiner.extract_bioactive_peptides(mbpdb)
        self.assertIsNotNone(grouped)
        # RELEELNVPGEI appears twice in mbpdb – should be grouped to one row
        self.assertEqual(len(grouped[grouped['search_peptide'] == 'RELEELNVPGEI']), 1)

    def test_functions_joined_with_semicolon(self):
        mbpdb = make_mbpdb_df()
        _, grouped = data_combiner.extract_bioactive_peptides(mbpdb)
        releelnvpgei_row = grouped[grouped['search_peptide'] == 'RELEELNVPGEI'].iloc[0]
        # 'Opioid' appears twice – after dedup it should be just 'Opioid'
        self.assertIn('Opioid', str(releelnvpgei_row['function']))

    def test_empty_mbpdb_returns_none(self):
        cleaned, grouped = data_combiner.extract_bioactive_peptides(None)
        self.assertIsNone(cleaned)
        self.assertIsNone(grouped)

    def test_none_protein_id_filtered(self):
        mbpdb = make_mbpdb_df().copy()
        mbpdb.loc[0, 'protein_id'] = 'None'
        cleaned, _ = data_combiner.extract_bioactive_peptides(mbpdb)
        self.assertNotIn('None', cleaned['protein_id'].values)


class TestCreateUniqueId(unittest.TestCase):

    def test_plain_sequence(self):
        row = pd.Series({'Sequence': 'RELEELNVPGEI', 'Modifications': None})
        uid = data_combiner.create_unique_id(row)
        self.assertEqual(uid, 'RELEELNVPGEI')

    def test_with_modifications(self):
        row = pd.Series({'Sequence': 'SEESITRINK', 'Modifications': '1xPhospho [S1]'})
        uid = data_combiner.create_unique_id(row)
        self.assertEqual(uid, 'SEESITRINK_1xPhospho [S1]')

    def test_no_trailing_underscore(self):
        row = pd.Series({'Sequence': 'ABCDE', 'Modifications': None})
        uid = data_combiner.create_unique_id(row)
        self.assertFalse(uid.endswith('_'))


class TestProcessPdResults(unittest.TestCase):

    def _get_validated_df(self):
        """Return a DataFrame as it would look after validation (Protein column renamed)."""
        df = make_peptidomic_df().rename(
            columns={'Master Protein Accessions': 'Protein'}
        )
        return df

    def test_adds_unique_peptide_id(self):
        df = self._get_validated_df()
        _, mbpdb_grouped = data_combiner.extract_bioactive_peptides(make_mbpdb_df())
        result = data_combiner.process_pd_results(df, mbpdb_grouped, make_protein_dict())
        self.assertIn('Unique Peptide ID', result.columns)

    def test_first_columns_order(self):
        """Unique Peptide ID, Sequence, Protein, Positions in Proteins come first."""
        df = self._get_validated_df()
        _, mbpdb_grouped = data_combiner.extract_bioactive_peptides(make_mbpdb_df())
        result = data_combiner.process_pd_results(df, mbpdb_grouped, make_protein_dict())
        expected_start = ['Unique Peptide ID', 'Sequence', 'Protein', 'Positions in Proteins']
        self.assertEqual(list(result.columns[:4]), expected_start)

    def test_merges_mbpdb_function(self):
        df = self._get_validated_df()
        _, mbpdb_grouped = data_combiner.extract_bioactive_peptides(make_mbpdb_df())
        result = data_combiner.process_pd_results(df, mbpdb_grouped, make_protein_dict())
        self.assertIn('function', result.columns)
        # RELEELNVPGEI should have 'Opioid' function
        row = result[result['Sequence'] == 'RELEELNVPGEI'].iloc[0]
        self.assertIn('Opioid', str(row['function']))

    def test_no_mbpdb_adds_nan_function(self):
        df = self._get_validated_df()
        result = data_combiner.process_pd_results(df, None, make_protein_dict())
        self.assertIn('function', result.columns)
        self.assertTrue(result['function'].isna().all())

    def test_start_end_positions_extracted(self):
        df = self._get_validated_df()
        result = data_combiner.process_pd_results(df, None, make_protein_dict())
        # P02666 [16-27] → start=16, end=27
        row = result[result['Sequence'] == 'RELEELNVPGEI'].iloc[0]
        self.assertEqual(row['start'], 16)
        self.assertEqual(row['end'], 27)


class TestCollapseMultiproteinDuplicateRows(unittest.TestCase):
    """Skyline reports a peptide mapped to multiple proteins as separate
    rows with otherwise-identical data (abundances included), unlike
    Proteome Discoverer/MaxQuant/PEAKS/Spectronaut which join multi-protein
    peptides into one row. Left uncollapsed this duplicates the same
    measured abundance under one 'Unique Peptide ID' — a real risk for any
    sum-based downstream analysis."""

    def test_identical_rows_collapsed_with_proteins_joined(self):
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK', 'PEPTIDEK'],
            'Protein': ['P12931', 'P07947'],
            'Protein Name': ['sp|P12931|SRC_HUMAN', 'sp|P07947|YES_HUMAN'],
            'Sample_A1': [100.0, 100.0],
        })
        result = data_combiner.process_pd_results(df, None, {})
        self.assertEqual(len(result), 1)
        row = result.iloc[0]
        self.assertEqual(row['Protein'], 'P12931; P07947')
        self.assertEqual(row['Protein Name'], 'sp|P12931|SRC_HUMAN; sp|P07947|YES_HUMAN')
        self.assertEqual(row['Sample_A1'], 100.0)

    def test_genuinely_different_abundance_not_collapsed(self):
        """If the 'duplicate' rows actually differ in a non-protein column,
        collapsing would silently discard real data — must be left as-is."""
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK', 'PEPTIDEK'],
            'Protein': ['P12931', 'P07947'],
            'Sample_A1': [100.0, 250.0],
        })
        result = data_combiner.process_pd_results(df, None, {})
        self.assertEqual(len(result), 2)

    def test_no_duplicates_unaffected(self):
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK', 'SAMPLER'],
            'Protein': ['P12931', 'P07947'],
            'Sample_A1': [100.0, 200.0],
        })
        result = data_combiner.process_pd_results(df, None, {})
        self.assertEqual(len(result), 2)

    def test_three_way_duplicate_all_proteins_joined(self):
        df = pd.DataFrame({
            'Sequence': ['PEPTIDEK'] * 3,
            'Protein': ['P1', 'P2', 'P3'],
            'Sample_A1': [50.0, 50.0, 50.0],
        })
        result = data_combiner.process_pd_results(df, None, {})
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]['Protein'], 'P1; P2; P3')


class TestCalculateGroupAbundanceAverages(unittest.TestCase):

    def test_avg_columns_added(self):
        df = make_peptidomic_df()
        gd = make_group_data()
        result = data_combiner.calculate_group_abundance_averages(df, gd)
        self.assertIn('Avg_GroupA', result.columns)
        self.assertIn('Avg_GroupB', result.columns)

    def test_avg_values_correct(self):
        df = make_peptidomic_df()
        gd = make_group_data()
        result = data_combiner.calculate_group_abundance_averages(df, gd)
        expected_a1 = (df['Sample_A1'] + df['Sample_A2']).values / 2
        np.testing.assert_allclose(result['Avg_GroupA'].values, expected_a1)

    def test_idempotent(self):
        """Calling twice should not duplicate columns."""
        df = make_peptidomic_df()
        gd = make_group_data()
        result = data_combiner.calculate_group_abundance_averages(df, gd)
        result2 = data_combiner.calculate_group_abundance_averages(result, gd)
        self.assertEqual(result2.columns.tolist(), result.columns.tolist())


class TestUpdateColumnNamesWithGroups(unittest.TestCase):

    def test_grouped_suffix_added(self):
        df = make_peptidomic_df()
        gd = make_group_data()
        result = data_combiner.update_column_names_with_groups(df, gd)
        # Sample_A1 belongs to GroupA
        self.assertIn("Sample_A1 'Grouped: (GroupA)'", result.columns)
        self.assertIn("Sample_B2 'Grouped: (GroupB)'", result.columns)


class TestProcessData(unittest.TestCase):
    """Integration test for the full data processing pipeline."""

    def test_full_pipeline_returns_dataframe(self):
        pd_results = make_peptidomic_df().rename(
            columns={'Master Protein Accessions': 'Protein'}
        )
        mbpdb = make_mbpdb_df()
        gd = make_group_data()
        protein_dict = make_protein_dict()

        result = data_combiner.process_data(pd_results, None, mbpdb, gd, protein_dict)

        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertGreater(len(result), 0)

    def test_full_pipeline_has_expected_columns(self):
        pd_results = make_peptidomic_df().rename(
            columns={'Master Protein Accessions': 'Protein'}
        )
        result = data_combiner.process_data(
            pd_results, None, make_mbpdb_df(), make_group_data(), make_protein_dict()
        )
        for col in ['Unique Peptide ID', 'Sequence', 'Protein',
                    'protein_species', 'protein_name', 'Avg_GroupA', 'Avg_GroupB']:
            self.assertIn(col, result.columns, msg=f"Missing column: {col}")

    def test_no_pd_results_returns_none(self):
        result = data_combiner.process_data(None, None, None, {}, {})
        self.assertIsNone(result)

    def test_empty_pd_results_returns_none(self):
        result = data_combiner.process_data(pd.DataFrame(), None, None, {}, {})
        self.assertIsNone(result)

    def test_without_mbpdb_still_works(self):
        pd_results = make_peptidomic_df().rename(
            columns={'Master Protein Accessions': 'Protein'}
        )
        result = data_combiner.process_data(
            pd_results, None, None, make_group_data(), make_protein_dict()
        )
        self.assertIsNotNone(result)

    def test_without_groups_still_works(self):
        pd_results = make_peptidomic_df().rename(
            columns={'Master Protein Accessions': 'Protein'}
        )
        result = data_combiner.process_data(
            pd_results, None, make_mbpdb_df(), None, make_protein_dict()
        )
        self.assertIsNotNone(result)


# ---------------------------------------------------------------------------
# export_manager tests
# ---------------------------------------------------------------------------

def _make_merged_df():
    """Build a merged DataFrame as process_data would produce it."""
    pd_results = make_peptidomic_df().rename(
        columns={'Master Protein Accessions': 'Protein'}
    )
    return data_combiner.process_data(
        pd_results, None, make_mbpdb_df(), make_group_data(), make_protein_dict()
    )


class TestExportMbpdbResults(unittest.TestCase):

    def test_returns_tsv_bytes(self):
        mbpdb = make_mbpdb_df()
        result = export_manager.export_mbpdb_results(mbpdb)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, bytes)
        # TSV should have tabs
        self.assertIn(b'\t', result)

    def test_none_returns_none(self):
        result = export_manager.export_mbpdb_results(None)
        self.assertIsNone(result)


class TestExportGroupDefinitions(unittest.TestCase):

    def test_returns_json_bytes(self):
        gd = make_group_data()
        result = export_manager.export_group_definitions(gd)
        self.assertIsNotNone(result)
        parsed = json.loads(result.decode())
        self.assertIn('GroupA', parsed)
        self.assertEqual(parsed['GroupA'], ['Sample_A1', 'Sample_A2'])

    def test_none_returns_none(self):
        self.assertIsNone(export_manager.export_group_definitions(None))
        self.assertIsNone(export_manager.export_group_definitions({}))


class TestExportMergedDataset(unittest.TestCase):

    def test_returns_csv_bytes(self):
        merged = _make_merged_df()
        result = export_manager.export_merged_dataset(merged)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, bytes)
        # CSV should have the Unique Peptide ID column header
        self.assertIn(b'Unique Peptide ID', result)

    def test_none_returns_none(self):
        self.assertIsNone(export_manager.export_merged_dataset(None))


class TestExportSequenceList(unittest.TestCase):

    def test_returns_csv_bytes(self):
        merged = _make_merged_df()
        result = export_manager.export_sequence_list(merged, make_group_data())
        # May be None if no non-zero averages; just check it doesn't crash
        if result is not None:
            self.assertIsInstance(result, bytes)


class TestSummedPeptideResults(unittest.TestCase):

    def test_returns_dict_with_groups(self):
        merged = _make_merged_df()
        result = export_manager.summed_peptide_results(merged, make_group_data())
        self.assertIsInstance(result, dict)
        self.assertIn('GroupA', result)
        self.assertIn('GroupB', result)

    def test_group_has_expected_keys(self):
        merged = _make_merged_df()
        result = export_manager.summed_peptide_results(merged, make_group_data())
        group_a = result['GroupA']
        for key in ['unique_peptides', 'total_Abundance', 'total_sem',
                    'abundance_sem', 'count_sem', 'replicate_data']:
            self.assertIn(key, group_a, msg=f"Missing key: {key}")

    def test_total_abundance_positive(self):
        merged = _make_merged_df()
        result = export_manager.summed_peptide_results(merged, make_group_data())
        self.assertGreater(result['GroupA']['total_Abundance'], 0)

    def test_export_summed_peptide_data_returns_xlsx(self):
        merged = _make_merged_df()
        result = export_manager.export_summed_peptide_data(merged, make_group_data())
        if result is not None:
            # XLSX files start with PK (zip magic bytes)
            self.assertEqual(result[:2], b'PK')


class TestExportProteinData(unittest.TestCase):

    def test_returns_xlsx_bytes(self):
        merged = _make_merged_df()
        result, warnings = export_manager.export_protein_data(
            merged, make_group_data(), make_protein_dict()
        )
        if result is not None:
            self.assertIsInstance(result, bytes)
            # XLSX files start with PK (zip magic bytes)
            self.assertEqual(result[:2], b'PK')

    def test_none_merged_returns_none(self):
        result, _ = export_manager.export_protein_data(None, make_group_data(), {})
        self.assertIsNone(result)


class TestExportGroupCorrelation(unittest.TestCase):

    def test_returns_xlsx(self):
        merged = _make_merged_df()
        result = export_manager.export_group_correlation(
            merged, make_group_data(), 'Pearson', True
        )
        self.assertIsNotNone(result)
        self.assertEqual(result[:2], b'PK')

    def test_spearman_works(self):
        merged = _make_merged_df()
        result = export_manager.export_group_correlation(
            merged, make_group_data(), 'Spearman', True
        )
        self.assertIsNotNone(result)


class TestExportReplicateCorrelation(unittest.TestCase):

    def test_returns_xlsx(self):
        merged = _make_merged_df()
        result = export_manager.export_replicate_correlation(
            merged, make_group_data(), 'Pearson', True
        )
        self.assertIsNotNone(result)
        self.assertEqual(result[:2], b'PK')


class TestCalculateCorrelation(unittest.TestCase):

    def test_pearson_perfect_correlation(self):
        x = pd.Series([1.0, 2.0, 3.0, 4.0])
        y = pd.Series([2.0, 4.0, 6.0, 8.0])
        corr, p_str = export_manager.calculate_correlation(x, y, 'Pearson')
        self.assertAlmostEqual(corr, 1.0, places=5)

    def test_spearman(self):
        x = pd.Series([1.0, 2.0, 3.0, 4.0])
        y = pd.Series([1.0, 4.0, 9.0, 16.0])
        corr, p_str = export_manager.calculate_correlation(x, y, 'Spearman')
        self.assertAlmostEqual(corr, 1.0, places=5)

    def test_p_value_format(self):
        x = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])
        y = pd.Series([1.1, 2.1, 3.1, 4.1, 5.1])
        _, p_str = export_manager.calculate_correlation(x, y)
        self.assertTrue(p_str.startswith('p'))


class TestPrepareData(unittest.TestCase):

    def test_log10_transform(self):
        data = pd.Series([10.0, 100.0, 1000.0])
        result = export_manager.prepare_data(data, log_transform=True)
        np.testing.assert_allclose(result.values, [1.0, 2.0, 3.0])

    def test_no_transform(self):
        data = pd.Series([10.0, 100.0])
        result = export_manager.prepare_data(data, log_transform=False)
        np.testing.assert_array_equal(result.values, data.values)

    def test_log2_transform(self):
        data = pd.Series([2.0, 8.0, 1024.0])
        result = export_manager.prepare_data(data, log_transform=True, log_base=2)
        np.testing.assert_allclose(result.values, [1.0, 3.0, 10.0])

    def test_log_label(self):
        self.assertEqual(export_manager.log_label(True, 2), 'Log2')
        self.assertEqual(export_manager.log_label(True, 10), 'Log10')
        self.assertEqual(export_manager.log_label(False), 'None')


# ---------------------------------------------------------------------------
# Notebook parity checks
# ---------------------------------------------------------------------------

class TestNotebookParity(unittest.TestCase):
    """
    Verify that the Django services produce the same logical output
    as the equivalent notebook class methods.
    """

    def test_unique_id_matches_notebook_logic(self):
        """
        Notebook create_unique_id: sequence + '_' + modifications (stripped),
        or just sequence if no modification. Trailing underscores removed.
        Django version should match.
        """
        cases = [
            # (sequence, modification, expected_id)
            ('RELEELNVPGEI', None, 'RELEELNVPGEI'),
            ('SEESITRINK', '1xPhospho [S1]', 'SEESITRINK_1xPhospho [S1]'),
            ('PEPTIDE', '', 'PEPTIDE'),
        ]
        for seq, mod, expected in cases:
            row = pd.Series({'Sequence': seq, 'Modifications': mod if mod else None})
            result = data_combiner.create_unique_id(row)
            self.assertEqual(result, expected, msg=f"Failed for seq={seq}, mod={mod}")

    def test_group_naming_format_matches_notebook(self):
        """
        Notebook update_column_names_with_groups renames columns to:
        "{column} 'Grouped: ({group_names})'"
        Django must match this format exactly.
        """
        df = pd.DataFrame({
            'Sample_A1': [1.0, 2.0],
            'Sample_A2': [3.0, 4.0],
        })
        gd = {'1': {'grouping_variable': 'GroupA', 'abundance_columns': ['Sample_A1', 'Sample_A2']}}
        result = data_combiner.update_column_names_with_groups(df, gd)
        self.assertIn("Sample_A1 'Grouped: (GroupA)'", result.columns)
        self.assertIn("Sample_A2 'Grouped: (GroupA)'", result.columns)

    def test_avg_column_naming_matches_notebook(self):
        """Average column name must be 'Avg_{grouping_variable}'."""
        df = pd.DataFrame({
            'Sample_A1': [1.0, 2.0],
            'Sample_A2': [3.0, 4.0],
        })
        gd = {'1': {'grouping_variable': 'MyGroup', 'abundance_columns': ['Sample_A1', 'Sample_A2']}}
        result = data_combiner.calculate_group_abundance_averages(df, gd)
        self.assertIn('Avg_MyGroup', result.columns)

    def test_mbpdb_function_aggregation_matches_notebook(self):
        """
        Notebook aggregates function as list of unique values joined with '; '.
        Django should produce the same output.
        """
        mbpdb = pd.DataFrame({
            'search_peptide': ['PEPTIDE', 'PEPTIDE'],
            'peptide':        ['PEPTIDE', 'PEPTIDE'],
            'protein_id':     ['P001', 'P001'],
            'function':       ['Opioid', 'ACE inhibitor'],
        })
        _, grouped = data_combiner.extract_bioactive_peptides(mbpdb)
        row = grouped[grouped['search_peptide'] == 'PEPTIDE'].iloc[0]
        functions = str(row['function'])
        self.assertIn('Opioid', functions)
        self.assertIn('ACE inhibitor', functions)
        # Notebook joins with '; '
        self.assertIn('; ', functions)

    def test_protein_info_column_position(self):
        """
        Notebook inserts protein_species and protein_name immediately after
        the protein accession column. Django must match this.
        """
        df = pd.DataFrame({
            'Unique Peptide ID': ['ABC'],
            'Sequence':          ['ABC'],
            'Protein':           ['P02666'],
            'Positions in Proteins': ['P02666 [1-3]'],
            'Sample_A1': [100.0],
        })
        result = data_combiner.add_protein_info(df, make_protein_dict())
        cols = list(result.columns)
        p_idx = cols.index('Protein')
        self.assertEqual(cols[p_idx + 1], 'protein_species')
        self.assertEqual(cols[p_idx + 2], 'protein_name')


# ---------------------------------------------------------------------------
# Real archive data tests — validate pipeline against actual lab files
# ---------------------------------------------------------------------------

ARCHIVE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    'examples', 'legacy_v1', 'archive'
)


def _load_real(filename, file_type='Peptidomic'):
    path = os.path.join(ARCHIVE, filename)
    with open(path, 'rb') as f:
        content = f.read()
    df, status, err, warn = data_loader.load_and_validate_file(
        content, filename, file_type
    )
    return df, status, err, warn


# ---------------------------------------------------------------------------
# not_in_use_examples — one real, unmodified export per engine (MaxQuant,
# Skyline, PEAKS, Spectronaut, plus a Skyline-based lab pipeline export)
# exercising Table 1's full column-mapping surface, including the
# long/per-run -> wide pivot (PEAKS, Spectronaut) and the Skyline
# multi-protein-row collapse. Tracked in git but not linked from any
# template — see the README in that directory.
# ---------------------------------------------------------------------------

ADDITIONAL_FORMATS = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    'examples', 'legacy_v1', 'additional_peptidomic_input_formats'
)


def _load_additional_format(filename, file_type='Peptidomic'):
    path = os.path.join(ADDITIONAL_FORMATS, filename)
    with open(path, 'rb') as f:
        content = f.read()
    df, status, err, warn = data_loader.load_and_validate_file(
        content, filename, file_type
    )
    return df, status, err, warn


def _run_full_pipeline(df):
    """Load -> detect abundance columns -> group (one group per column) ->
    process. Mirrors the wizard's Step 1 -> Step 2 -> Step 4 path so a real
    export's *output* shape can be checked, not just that it loads."""
    filtered = data_loader.get_filtered_columns(df)
    if len(filtered) < 2:
        filtered = data_loader.get_filtered_columns_fallback(df)
    gd, err = group_processing.build_no_group_data(filtered)
    result = data_combiner.process_data(df, None, None, gd, {})
    return result, filtered


class TestRealData_MaxQuant(unittest.TestCase):
    """max_quant.csv — real MaxQuant peptides.txt export, wide (one
    'Intensity <expt>' column per sample), 25 rows, no pivot needed."""

    FILE = 'max_quant.csv'

    def setUp(self):
        self.df, self.status, self.err, self.warn = _load_additional_format(self.FILE)

    def test_loads_successfully(self):
        self.assertEqual(self.status, 'yes', msg=self.err)

    def test_not_treated_as_long_format(self):
        self.assertNotIn('long-format', self.warn or '')

    def test_protein_and_sequence_mapped(self):
        self.assertIn('Protein', self.df.columns)
        self.assertIn('Sequence', self.df.columns)

    def test_full_pipeline_matches_reference_schema(self):
        result, filtered = _run_full_pipeline(self.df)
        self.assertIsNotNone(result)
        self.assertEqual(len(result), len(self.df))  # already peptide-per-row
        for col in ('Unique Peptide ID', 'Sequence', 'Protein', 'protein_species', 'protein_name'):
            self.assertIn(col, result.columns)
        self.assertEqual(result['Unique Peptide ID'].duplicated().sum(), 0)
        self.assertTrue(any(c.startswith('Avg_') for c in result.columns))


class TestRealData_SkylineExample(unittest.TestCase):
    """'Skyline example.csv' — real Skyline custom report, wide (one
    'Normalized Area' column per sample). Reports a peptide mapped to
    multiple proteins as separate rows with identical data, which
    process_pd_results must collapse to one row per peptide."""

    FILE = 'Skyline example.csv'

    def setUp(self):
        self.df, self.status, self.err, self.warn = _load_additional_format(self.FILE)

    def test_loads_successfully(self):
        self.assertEqual(self.status, 'yes', msg=self.err)

    def test_not_treated_as_long_format(self):
        self.assertNotIn('long-format', self.warn or '')

    def test_raw_file_has_multiprotein_duplicate_rows(self):
        """Sanity check on the fixture itself: confirms the scenario this
        test class exists to cover is actually present in the real file."""
        self.assertGreater(len(self.df), self.df['Peptide Modified Sequence'].nunique())

    def test_full_pipeline_collapses_multiprotein_rows(self):
        result, filtered = _run_full_pipeline(self.df)
        self.assertIsNotNone(result)
        self.assertGreater(len(filtered), 0)
        # One row per Unique Peptide ID — no duplicates left over from
        # Skyline's one-row-per-(peptide, protein) convention.
        self.assertEqual(result['Unique Peptide ID'].duplicated().sum(), 0)
        for col in ('Unique Peptide ID', 'Sequence', 'Protein', 'protein_species', 'protein_name'):
            self.assertIn(col, result.columns)
        # A peptide known (from the raw file) to map to 2 proteins should
        # now appear once, with both proteins joined.
        row = result[result['Unique Peptide ID'] == 'AANILVGENLVC[+57]K']
        self.assertEqual(len(row), 1)
        self.assertIn(';', str(row.iloc[0]['Protein']))


class TestRealData_PepexInput(unittest.TestCase):
    """PEPEX_input.tsv — the lab's own Skyline-based pipeline export, wide
    (one 'Total Area' column per sample), 993 rows, no pivot needed."""

    FILE = 'PEPEX_input.tsv'

    def setUp(self):
        self.df, self.status, self.err, self.warn = _load_additional_format(self.FILE)

    def test_loads_successfully(self):
        self.assertEqual(self.status, 'yes', msg=self.err)

    def test_not_treated_as_long_format(self):
        self.assertNotIn('long-format', self.warn or '')

    def test_full_pipeline_matches_reference_schema(self):
        result, filtered = _run_full_pipeline(self.df)
        self.assertIsNotNone(result)
        self.assertEqual(len(result), len(self.df))
        for col in ('Unique Peptide ID', 'Sequence', 'Protein', 'protein_species', 'protein_name'):
            self.assertIn(col, result.columns)
        self.assertEqual(result['Unique Peptide ID'].duplicated().sum(), 0)


class TestRealData_PeaksLongFormat(unittest.TestCase):
    """PEAKS_example.csv — real PEAKS peptide.csv export (immunopeptidomics/
    MAPPs style, timsTOF ion-mobility columns), one row per PSM per
    'Source File'. Must be auto-pivoted before the rest of the pipeline
    (moved here from needs_transformation/ now that pivoting is automatic)."""

    FILE = 'PEAKS_example.csv'

    def setUp(self):
        self.df, self.status, self.err, self.warn = _load_additional_format(self.FILE)

    def test_loads_successfully(self):
        self.assertEqual(self.status, 'yes', msg=self.err)

    def test_pivoted_to_one_row_per_peptide(self):
        self.assertIn('long-format', self.warn or '')
        self.assertEqual(len(self.df), self.df['Sequence'].nunique())

    def test_full_pipeline_matches_reference_schema(self):
        result, filtered = _run_full_pipeline(self.df)
        self.assertIsNotNone(result)
        self.assertEqual(len(filtered), 8)  # 8 raw files (samples)
        for col in ('Unique Peptide ID', 'Sequence', 'Protein', 'protein_species', 'protein_name'):
            self.assertIn(col, result.columns)
        self.assertEqual(result['Unique Peptide ID'].duplicated().sum(), 0)
        self.assertTrue(any(c.startswith('Avg_') for c in result.columns))


class TestRealData_SpectronautLongFormat(unittest.TestCase):
    """spectronaut.tsv — real Spectronaut Normal Report (DIA, HYE mix
    conditions), one row per precursor per 'R.FileName'. Must be
    auto-pivoted (moved here from needs_transformation/ now that pivoting
    is automatic)."""

    FILE = 'spectronaut.tsv'

    def setUp(self):
        self.df, self.status, self.err, self.warn = _load_additional_format(self.FILE)

    def test_loads_successfully(self):
        self.assertEqual(self.status, 'yes', msg=self.err)

    def test_pivoted_to_one_row_per_peptide(self):
        self.assertIn('long-format', self.warn or '')

    def test_full_pipeline_matches_reference_schema(self):
        result, filtered = _run_full_pipeline(self.df)
        self.assertIsNotNone(result)
        self.assertEqual(len(filtered), 8)  # 8 raw files (samples)
        for col in ('Unique Peptide ID', 'Sequence', 'Protein', 'protein_species', 'protein_name'):
            self.assertIn(col, result.columns)
        self.assertEqual(result['Unique Peptide ID'].duplicated().sum(), 0)
        self.assertTrue(any(c.startswith('Avg_') for c in result.columns))


class TestRealData_ProteomeDiscoverer(unittest.TestCase):
    """test.csv — Proteome Discoverer format (Annotated Sequence, Master Protein Accessions)"""

    FILE = 'test.csv'

    def setUp(self):
        self.df, self.status, self.err, self.warn = _load_real(self.FILE)

    def test_loads_successfully(self):
        self.assertEqual(self.status, 'yes', msg=self.err)
        self.assertIsNotNone(self.df)

    def test_protein_column_mapped(self):
        self.assertIn('Protein', self.df.columns)

    def test_sequence_present(self):
        self.assertTrue(
            'Sequence' in self.df.columns or 'Annotated Sequence' in self.df.columns
        )

    def test_has_rows(self):
        self.assertGreater(len(self.df), 0)

    def test_filtered_columns_are_abundance(self):
        filtered = data_loader.get_filtered_columns(self.df)
        self.assertGreater(len(filtered), 0)
        # None of the filtered columns should be metadata
        for col in filtered:
            self.assertNotIn(col, ['Sequence', 'Protein', 'Modifications',
                                   'Annotated Sequence', 'Positions in Proteins'])

    def test_full_pipeline(self):
        gd = {}
        filtered = data_loader.get_filtered_columns(self.df)
        if filtered:
            gd, _ = group_processing.build_no_group_data(filtered)
        result = data_combiner.process_data(self.df, None, None, gd, {})
        self.assertIsNotNone(result)
        self.assertIn('Unique Peptide ID', result.columns)
        self.assertIn('Sequence', result.columns)


class TestRealData_ExampleShort(unittest.TestCase):
    """example_peptide_data_short.csv — PD format with direct Sequence column, ~49 rows"""

    FILE = 'example_peptide_data_short.csv'

    def setUp(self):
        self.df, self.status, self.err, self.warn = _load_real(self.FILE)

    def test_loads_successfully(self):
        self.assertEqual(self.status, 'yes', msg=self.err)

    def test_protein_column_mapped(self):
        self.assertIn('Protein', self.df.columns)

    def test_abundance_columns_detected(self):
        filtered = data_loader.get_filtered_columns(self.df)
        self.assertGreater(len(filtered), 0)

    def test_sequences_extractable(self):
        seqs = data_loader.extract_sequences(self.df)
        self.assertGreater(len(seqs), 0)
        # All should be non-empty strings
        for s in seqs:
            self.assertIsInstance(s, str)
            self.assertGreater(len(s), 0)

    def test_full_pipeline_with_groups(self):
        filtered = data_loader.get_filtered_columns(self.df)
        self.assertGreater(len(filtered), 0)
        # Put all columns into one group
        gd, err = group_processing.build_group_data(filtered[:3], 'TestGroup')
        self.assertIsNone(err)
        result = data_combiner.process_data(self.df, None, None, gd, {})
        self.assertIsNotNone(result)
        self.assertIn('Avg_TestGroup', result.columns)

    def test_export_merged_dataset(self):
        filtered = data_loader.get_filtered_columns(self.df)
        gd, _ = group_processing.build_no_group_data(filtered[:2])
        result = data_combiner.process_data(self.df, None, None, gd, {})
        csv_bytes = export_manager.export_merged_dataset(result)
        self.assertIsNotNone(csv_bytes)
        # Parse the CSV back and verify shape
        result_back = pd.read_csv(io.BytesIO(csv_bytes))
        self.assertEqual(len(result_back), len(result))


class TestRealData_Skyline(unittest.TestCase):
    """skyline_raw.csv — Skyline format (Protein, Peptide Modified Sequence, Peptide start/end), ~993 rows"""

    FILE = 'skyline_raw.csv'

    def setUp(self):
        self.df, self.status, self.err, self.warn = _load_real(self.FILE)

    def test_loads_successfully(self):
        self.assertEqual(self.status, 'yes', msg=self.err)

    def test_protein_column_present(self):
        self.assertIn('Protein', self.df.columns)

    def test_abundance_columns_detected(self):
        filtered = data_loader.get_filtered_columns(self.df)
        self.assertGreater(len(filtered), 0)

    def test_sequences_extractable(self):
        seqs = data_loader.extract_sequences(self.df)
        self.assertGreater(len(seqs), 0)

    def test_unique_id_no_trailing_underscore(self):
        """Every unique ID produced from real Skyline data must not end with '_'."""
        filtered = data_loader.get_filtered_columns(self.df)
        gd, _ = group_processing.build_no_group_data(filtered[:1])
        result = data_combiner.process_data(self.df, None, None, gd, {})
        if result is not None:
            for uid in result['Unique Peptide ID'].dropna():
                self.assertFalse(str(uid).endswith('_'), msg=f"Trailing underscore: {uid!r}")

    def test_full_pipeline_large_file(self):
        filtered = data_loader.get_filtered_columns(self.df)
        self.assertGreater(len(filtered), 0)
        gd, _ = group_processing.build_no_group_data(filtered[:4])
        result = data_combiner.process_data(self.df, None, None, gd, {})
        self.assertIsNotNone(result)
        # Row count should be preserved
        self.assertGreater(len(result), 0)


class TestRealData_Noom(unittest.TestCase):
    """noom_peptide_data.csv — Annotated Sequence format, ~2366 rows, large real dataset"""

    FILE = 'noom_peptide_data.csv'

    def setUp(self):
        self.df, self.status, self.err, self.warn = _load_real(self.FILE)

    def test_loads_successfully(self):
        self.assertEqual(self.status, 'yes', msg=self.err)

    def test_abundance_columns_detected(self):
        filtered = data_loader.get_filtered_columns(self.df)
        self.assertGreater(len(filtered), 4,
            msg=f"Expected multiple abundance cols, got: {filtered}")

    def test_group_json_round_trip(self):
        """Build groups from first 4 abundance columns, export JSON, re-import."""
        filtered = data_loader.get_filtered_columns(self.df)
        cols_a = filtered[:2]
        cols_b = filtered[2:4]
        gd = {}
        gd, _ = group_processing.build_group_data(cols_a, 'GroupA', gd)
        gd, _ = group_processing.build_group_data(cols_b, 'GroupB', gd)

        json_str = group_processing.export_group_data_json(gd)
        gd2, err = group_processing.parse_group_json(
            json_str.encode(), data_loader.get_filtered_columns(self.df)
        )
        self.assertIsNone(err)
        self.assertEqual(
            gd['GroupA']['abundance_columns'],
            gd2['GroupA']['abundance_columns']
        )

    def test_full_pipeline_produces_averages(self):
        filtered = data_loader.get_filtered_columns(self.df)
        gd, _ = group_processing.build_group_data(filtered[:3], 'GroupA')
        gd, _ = group_processing.build_group_data(filtered[3:6], 'GroupB', gd)
        result = data_combiner.process_data(self.df, None, None, gd, {})
        self.assertIsNotNone(result)
        self.assertIn('Avg_GroupA', result.columns)
        self.assertIn('Avg_GroupB', result.columns)
        # Averages should be numeric
        self.assertTrue(pd.api.types.is_numeric_dtype(result['Avg_GroupA']))

    def test_summed_peptide_export(self):
        filtered = data_loader.get_filtered_columns(self.df)
        gd, _ = group_processing.build_group_data(filtered[:3], 'GroupA')
        gd, _ = group_processing.build_group_data(filtered[3:6], 'GroupB', gd)
        result = data_combiner.process_data(self.df, None, None, gd, {})
        xlsx = export_manager.export_summed_peptide_data(result, gd)
        self.assertIsNotNone(xlsx)
        self.assertEqual(xlsx[:2], b'PK')  # XLSX is a zip

    def test_replicate_correlation_export(self):
        filtered = data_loader.get_filtered_columns(self.df)
        gd, _ = group_processing.build_group_data(filtered[:3], 'GroupA')
        result = data_combiner.process_data(self.df, None, None, gd, {})
        xlsx = export_manager.export_replicate_correlation(result, gd)
        self.assertIsNotNone(xlsx)
        self.assertEqual(xlsx[:2], b'PK')


class TestRealData_ManyPositions(unittest.TestCase):
    """peptides_with_many_positions.csv — rows with multiple protein positions ('; ' separated)"""

    FILE = 'peptides_with_many_positions.csv'

    def setUp(self):
        self.df, self.status, self.err, self.warn = _load_real(self.FILE)

    def test_loads_successfully(self):
        self.assertEqual(self.status, 'yes', msg=self.err)

    def test_multiprotein_rows_handled(self):
        """Rows with '; ' in Positions in Proteins should not crash processing."""
        if 'Positions in Proteins' in self.df.columns:
            multi = self.df['Positions in Proteins'].astype(str).str.contains('; ', na=False)
            self.assertGreater(multi.sum(), 0,
                msg="Expected rows with multiple protein positions")

        filtered = data_loader.get_filtered_columns(self.df)
        if not filtered:
            self.skipTest("No abundance columns found")
        gd, _ = group_processing.build_no_group_data(filtered[:2])
        result = data_combiner.process_data(self.df, None, None, gd, {})
        self.assertIsNotNone(result)

    def test_missing_proteins_detected(self):
        missing = data_loader.find_missing_proteins(self.df, {})
        # With empty protein_dict every protein should be "missing"
        if 'Protein' in self.df.columns:
            self.assertGreater(len(missing), 0)

    def test_protein_combinations_extracted(self):
        from data_transformation.services import protein_handler
        combos = protein_handler.get_protein_combinations(self.df, {})
        # File has multi-protein rows — should find at least one combination
        self.assertIsInstance(combos, list)


class TestDecisionsToUi(unittest.TestCase):
    """Restore-hints for the Step-3 protein-mapping UI.

    Regression coverage for the bug where a saved Custom ID (and split/asis
    selections) were not restored when returning to the Protein Mapping step.
    """

    def setUp(self):
        from data_transformation.views import _decisions_to_ui
        self._to_ui = _decisions_to_ui

    def test_ui_custom_id_restored(self):
        out = self._to_ui({'P02666A1; P02666A2': {'action': 'CUSTOM', 'protein_id': 'P02666'}})
        self.assertEqual(out['P02666A1; P02666A2'],
                         {'mode': 'custom', 'protein_id': 'P02666'})

    def test_ui_split_selection_restored(self):
        out = self._to_ui({'Q1; Q2; Q3': {'action': 'SPLIT', 'protein_ids': ['Q1', 'Q3']}})
        self.assertEqual(out['Q1; Q2; Q3'], {'mode': 'split', 'protein_ids': ['Q1', 'Q3']})

    def test_ui_asis_restored(self):
        out = self._to_ui({'R1; R2': {'action': 'ASIS'}})
        self.assertEqual(out['R1; R2'], {'mode': 'asis'})

    def test_backend_format_custom_restored(self):
        # Format produced by an uploaded mapping key
        out = self._to_ui({'S1; S2': {'S1': 'CUSTOM:MYID', 'S2': 'REMOVE'}})
        self.assertEqual(out['S1; S2'], {'mode': 'custom', 'protein_id': 'MYID'})

    def test_custom_name_restored(self):
        out = self._to_ui({'P02666A1; P02666A2': {
            'action': 'CUSTOM', 'protein_id': 'P02666', 'protein_name': 'Beta-casein'}})
        self.assertEqual(out['P02666A1; P02666A2'],
                         {'mode': 'custom', 'protein_id': 'P02666', 'protein_name': 'Beta-casein'})

    def test_custom_without_name_omits_key(self):
        # Unchanged behavior when no name was entered.
        out = self._to_ui({'P02666A1; P02666A2': {'action': 'CUSTOM', 'protein_id': 'P02666'}})
        self.assertNotIn('protein_name', out['P02666A1; P02666A2'])

    def test_backend_format_split_restored(self):
        out = self._to_ui({'T1; T2': {'T1': 'NEW', 'T2': 'NEW'}})
        self.assertEqual(out['T1; T2'], {'mode': 'split', 'protein_ids': ['T1', 'T2']})

    def test_empty_and_malformed_ignored(self):
        self.assertEqual(self._to_ui({}), {})
        self.assertEqual(self._to_ui({'X': 'not-a-dict'}), {})


class TestApplyCustomNames(unittest.TestCase):
    """Custom protein names entered alongside a CUSTOM decision are persisted
    into protein_dict so add_protein_info fills them in at process time."""

    def setUp(self):
        from data_transformation.views import _apply_custom_names
        self._apply = _apply_custom_names

    def test_custom_name_written_to_protein_dict(self):
        protein_dict = {}
        self._apply(
            {'P02666A1; P02666A2': {'action': 'CUSTOM', 'protein_id': 'P02666',
                                     'protein_name': 'Beta-casein'}},
            protein_dict)
        self.assertEqual(protein_dict['P02666']['name'], 'Beta-casein')

    def test_custom_without_name_leaves_dict_untouched(self):
        protein_dict = {}
        self._apply({'A; B': {'action': 'CUSTOM', 'protein_id': 'X'}}, protein_dict)
        self.assertEqual(protein_dict, {})

    def test_non_custom_actions_ignored(self):
        protein_dict = {}
        self._apply({'A; B': {'action': 'SPLIT', 'protein_ids': ['A']}}, protein_dict)
        self.assertEqual(protein_dict, {})

    def test_legacy_backend_format_ignored(self):
        # No 'action' key → legacy backend format, which has no room for a name.
        protein_dict = {}
        self._apply({'A; B': {'A': 'CUSTOM:X'}}, protein_dict)
        self.assertEqual(protein_dict, {})

    def test_existing_entry_updated_not_overwritten(self):
        protein_dict = {'P02666': {'species': 'Bovine'}}
        self._apply(
            {'A; B': {'action': 'CUSTOM', 'protein_id': 'P02666', 'protein_name': 'Beta-casein'}},
            protein_dict)
        self.assertEqual(protein_dict['P02666'], {'species': 'Bovine', 'name': 'Beta-casein'})


def make_beta_casein_df():
    """
    Mirrors the β-casein variant situation in merged_dataframetemp.csv:
    a multi-accession row already collapsed to P02666, plus singleton
    P02666A1 / P02666A2 rows that the combination step can't reach.
    """
    return pd.DataFrame({
        'Sequence': ['PPFLQPEVM', 'AAAA', 'BBBB', 'CCCC'],
        'Protein': ['P02666', 'P02666A2', 'P02666A1', 'P08037'],
        'protein_species': ['Bovine', 'Bovine', 'Bovine', 'Bovine'],
        'protein_name': [
            'CASB_BOVIN Beta-casein',
            'CASB_BOVIN Beta-casein A2',
            'CASB_BOVIN Beta-casein A1',
            'B4GT1_BOVIN Beta-1,4-galactosyltransferase 1',
        ],
        'Positions in Proteins': [
            'P02666 [85-93]',
            'P02666A2 [61-68]',
            'P02666A1 [61-69]',
            'P08037 [68-79]',
        ],
        'start': [85, 61, 61, 68],
        'end': [93, 68, 69, 79],
        'Master Protein Accessions': ['P02666A1; P02666A2', 'P02666A2', 'P02666A1', 'P08037'],
    })


class TestGetProteinSources(unittest.TestCase):
    def test_lists_every_accession_with_counts(self):
        df = make_beta_casein_df()
        protein_dict = {'P02666': {'name': 'Beta-casein', 'species': 'Bovine'}}
        sources = protein_handler.get_protein_sources(df, protein_dict)
        by_id = {s['id']: s for s in sources}
        self.assertEqual(set(by_id), {'P02666', 'P02666A1', 'P02666A2', 'P08037'})
        self.assertEqual(by_id['P02666']['count'], 1)
        self.assertEqual(by_id['P02666A1']['count'], 1)
        # dict-backed name/species come through for known accessions
        self.assertEqual(by_id['P02666']['name'], 'Beta-casein')

    def test_counts_multi_protein_members(self):
        df = pd.DataFrame({
            'Protein': ['P1; P2', 'P1'],
            'Positions in Proteins': ['P1 [1-5]; P2 [3-7]', 'P1 [1-5]'],
        })
        sources = protein_handler.get_protein_sources(df, {})
        by_id = {s['id']: s['count'] for s in sources}
        self.assertEqual(by_id, {'P1': 2, 'P2': 1})

    def test_empty_df(self):
        self.assertEqual(protein_handler.get_protein_sources(pd.DataFrame(), {}), [])


class TestApplySourceRenames(unittest.TestCase):
    def _merge_group(self):
        return [{
            'sources': ['P02666', 'P02666A1', 'P02666A2'],
            'target_id': 'P02666',
            'target_name': 'Beta-casein',
        }]

    def test_variants_unified_to_single_target(self):
        df = make_beta_casein_df()
        protein_dict = {}
        out, err = protein_handler.apply_source_renames(df, self._merge_group(), protein_dict)
        self.assertIsNone(err)
        beta = out[out['Sequence'].isin(['PPFLQPEVM', 'AAAA', 'BBBB'])]
        # All three β-casein rows collapse to the canonical accession + name
        self.assertTrue((beta['Protein'] == 'P02666').all())
        self.assertTrue((beta['protein_name'] == 'Beta-casein').all())
        # Species is not tracked by merge/rename groups; existing per-row
        # values are left untouched (all 'Bovine' in this fixture already).
        self.assertTrue((beta['protein_species'] == 'Bovine').all())
        # Positions rewritten to the target prefix, ranges (and start/end) retained
        pos = dict(zip(out['Sequence'], out['Positions in Proteins']))
        self.assertEqual(pos['AAAA'], 'P02666 [61-68]')
        self.assertEqual(pos['BBBB'], 'P02666 [61-69]')
        start = dict(zip(out['Sequence'], out['start']))
        self.assertEqual(start['AAAA'], 61)
        self.assertEqual(start['BBBB'], 61)

    def test_master_protein_accessions_untouched(self):
        df = make_beta_casein_df()
        out, err = protein_handler.apply_source_renames(df, self._merge_group(), {})
        self.assertIsNone(err)
        mpa = dict(zip(out['Sequence'], out['Master Protein Accessions']))
        self.assertEqual(mpa['PPFLQPEVM'], 'P02666A1; P02666A2')
        self.assertEqual(mpa['AAAA'], 'P02666A2')
        self.assertEqual(mpa['BBBB'], 'P02666A1')

    def test_unrelated_protein_unchanged(self):
        df = make_beta_casein_df()
        out, _ = protein_handler.apply_source_renames(df, self._merge_group(), {})
        row = out[out['Sequence'] == 'CCCC'].iloc[0]
        self.assertEqual(row['Protein'], 'P08037')
        self.assertEqual(row['Positions in Proteins'], 'P08037 [68-79]')

    def test_member_in_multi_protein_row_replaced_and_deduped(self):
        df = pd.DataFrame({
            'Sequence': ['X', 'Y'],
            'Protein': ['P02666A1; P02668', 'P02666A1; P02666A2'],
            'protein_name': ['n1', 'n2'],
            'protein_species': ['Bovine', 'Bovine'],
            'Positions in Proteins': [
                'P02666A1 [61-69]; P02668 [109-122]',
                'P02666A1 [61-69]; P02666A2 [61-68]',
            ],
        })
        groups = [{'sources': ['P02666A1', 'P02666A2'], 'target_id': 'P02666',
                   'target_name': 'Beta-casein'}]
        out, err = protein_handler.apply_source_renames(df, groups, {})
        self.assertIsNone(err)
        rows = dict(zip(out['Sequence'], out['Protein']))
        # In-place member replacement, still multi-protein → name/species left alone
        self.assertEqual(rows['X'], 'P02666; P02668')
        # Two variants collapse to one accession, de-duplicated to a single-protein row
        self.assertEqual(rows['Y'], 'P02666')
        yrow = out[out['Sequence'] == 'Y'].iloc[0]
        self.assertEqual(yrow['protein_name'], 'Beta-casein')
        # Positions de-duplicated where identical after rewrite is not the case here,
        # but the P02666A2 part rewrites to a distinct range so both remain
        self.assertIn('P02666 [61-69]', yrow['Positions in Proteins'])
        self.assertIn('P02666 [61-68]', yrow['Positions in Proteins'])

    def test_target_name_falls_back_to_protein_dict(self):
        df = make_beta_casein_df()
        protein_dict = {'P02666': {'name': 'Dict Beta-casein', 'species': 'Cow'}}
        groups = [{'sources': ['P02666A1'], 'target_id': 'P02666', 'target_name': ''}]
        out, err = protein_handler.apply_source_renames(df, groups, protein_dict)
        self.assertIsNone(err)
        row = out[out['Sequence'] == 'BBBB'].iloc[0]
        self.assertEqual(row['protein_name'], 'Dict Beta-casein')
        # Species is never written by merge/rename groups, even via protein_dict.
        self.assertEqual(row['protein_species'], 'Bovine')

    def test_species_column_never_touched(self):
        df = make_beta_casein_df()
        protein_dict = {'P02666': {'name': 'Beta-casein', 'species': 'SomeOtherSpecies'}}
        out, err = protein_handler.apply_source_renames(df, self._merge_group(), protein_dict)
        self.assertIsNone(err)
        beta = out[out['Sequence'].isin(['PPFLQPEVM', 'AAAA', 'BBBB'])]
        self.assertTrue((beta['protein_species'] == 'Bovine').all())

    def test_validation_errors(self):
        df = make_beta_casein_df()
        _, err = protein_handler.apply_source_renames(
            df, [{'sources': [], 'target_id': 'P02666'}], {})
        self.assertIsNotNone(err)
        _, err = protein_handler.apply_source_renames(
            df, [{'sources': ['P02666A1'], 'target_id': ''}], {})
        self.assertIsNotNone(err)
        # Same source pointed at two targets
        _, err = protein_handler.apply_source_renames(df, [
            {'sources': ['P02666A1'], 'target_id': 'P02666'},
            {'sources': ['P02666A1'], 'target_id': 'Q99999'},
        ], {})
        self.assertIsNotNone(err)


class TestDetectRenameConflicts(unittest.TestCase):
    def test_overlap_with_non_asis_decision_warns(self):
        renames = [{'sources': ['P02666A1', 'P02666A2'], 'target_id': 'P02666'}]
        decisions = {'P02666A1; P02666A2': {'action': 'CUSTOM', 'protein_id': 'P02666'}}
        warnings = protein_handler.detect_rename_conflicts(renames, decisions)
        self.assertEqual(len(warnings), 1)
        self.assertIn('P02666A1', warnings[0])

    def test_backend_format_overlap_warns(self):
        renames = [{'sources': ['S1'], 'target_id': 'S1'}]
        decisions = {'S1; S2': {'S1': 'CUSTOM:MYID', 'S2': 'REMOVE'}}
        warnings = protein_handler.detect_rename_conflicts(renames, decisions)
        self.assertEqual(len(warnings), 1)

    def test_asis_decision_does_not_warn(self):
        renames = [{'sources': ['P02666A1'], 'target_id': 'P02666'}]
        decisions = {'P02666A1; P02666A2': {'action': 'ASIS'}}
        self.assertEqual(protein_handler.detect_rename_conflicts(renames, decisions), [])

    def test_no_overlap_no_warning(self):
        renames = [{'sources': ['P02666A1'], 'target_id': 'P02666'}]
        decisions = {'Q1; Q2': {'action': 'SPLIT', 'protein_ids': ['Q1']}}
        self.assertEqual(protein_handler.detect_rename_conflicts(renames, decisions), [])


class TestMergeGroupsAsDecisions(unittest.TestCase):
    """Unified protein-mapping-key format: merge/rename groups folded into
    protein_decisions as action:"MERGE" entries, keyed like other combo
    decisions (sorted, semicolon-joined source accessions)."""

    def test_group_becomes_merge_entry(self):
        groups = [{'sources': ['P02666A2', 'P02666A1'], 'target_id': 'P02666',
                   'target_name': 'Beta-casein'}]
        entries = protein_handler.merge_groups_as_decisions(groups)
        self.assertEqual(entries, {
            'P02666A1; P02666A2': {
                'action': 'MERGE',
                'protein_ids': ['P02666A2', 'P02666A1'],
                'target_id': 'P02666',
                'target_name': 'Beta-casein',
            }
        })

    def test_empty_and_invalid_groups_produce_no_entries(self):
        self.assertEqual(protein_handler.merge_groups_as_decisions([]), {})
        self.assertEqual(protein_handler.merge_groups_as_decisions(None), {})
        # Invalid group (no target) is dropped rather than raising here —
        # validation errors surface through apply_source_renames instead.
        self.assertEqual(
            protein_handler.merge_groups_as_decisions([{'sources': ['A'], 'target_id': ''}]), {})

    def test_round_trip_through_extract(self):
        groups = [{'sources': ['P02666A1', 'P02666A2'], 'target_id': 'P02666',
                   'target_name': 'Beta-casein'}]
        entries = protein_handler.merge_groups_as_decisions(groups)
        combined = {'Q1; Q2': {'action': 'SPLIT', 'protein_ids': ['Q1']}, **entries}
        remaining, extracted_groups = protein_handler.extract_merge_groups(combined)
        self.assertEqual(remaining, {'Q1; Q2': {'action': 'SPLIT', 'protein_ids': ['Q1']}})
        self.assertEqual(extracted_groups, [
            {'sources': ['P02666A1', 'P02666A2'], 'target_id': 'P02666',
             'target_name': 'Beta-casein'}
        ])

    def test_extract_falls_back_to_key_when_protein_ids_missing(self):
        combined = {'P1; P2': {'action': 'MERGE', 'target_id': 'P1', 'target_name': ''}}
        remaining, groups = protein_handler.extract_merge_groups(combined)
        self.assertEqual(remaining, {})
        self.assertEqual(groups, [{'sources': ['P1', 'P2'], 'target_id': 'P1', 'target_name': ''}])

    def test_extract_leaves_non_merge_entries_untouched(self):
        combined = {
            'A; B': {'action': 'ASIS'},
            'legacy': {'X': 'CUSTOM:Y'},
        }
        remaining, groups = protein_handler.extract_merge_groups(combined)
        self.assertEqual(remaining, combined)
        self.assertEqual(groups, [])

    def test_empty_and_none_input(self):
        self.assertEqual(protein_handler.extract_merge_groups({}), ({}, []))
        self.assertEqual(protein_handler.extract_merge_groups(None), ({}, []))


if __name__ == '__main__':
    unittest.main()
