"""Ad hoc verification: column mapping is case-insensitive across example
formats, including deliberately case-mangled headers. Not part of the
permanent suite; scratch check for the Table 2 case-sensitivity fix."""
import io
import pandas as pd
import pytest

from data_transformation.services import data_loader

EXAMPLES = "peptide/static/peptide/examples/not_in_use_examples"

FILES = {
    "max_quant.csv": (",", "Sequence"),
    "PEAKS_example.csv": (",", "Peptide"),
    "Skyline example.csv": (",", "Peptide"),
    "spectronaut.tsv": ("\t", "PEP.StrippedSequence"),
    "PEPEX_input.tsv": ("\t", "Peptide"),
}


def _read(fname, sep):
    return pd.read_csv(f"{EXAMPLES}/{fname}", sep=sep, nrows=200)


def _to_bytes(df, sep):
    buf = io.BytesIO()
    df.to_csv(buf, index=False, sep=sep)
    return buf.getvalue()


@pytest.mark.parametrize("fname,sep,seq_col", [(f, s, c) for f, (s, c) in FILES.items()])
def test_original_casing(fname, sep, seq_col):
    df = _read(fname, sep)
    content = _to_bytes(df, sep)
    result_df, status, err, warn = data_loader.load_and_validate_file(content, fname, "Peptidomic")
    assert status == "yes", f"{fname}: {err}"
    assert "Sequence" in result_df.columns
    assert "Protein" in result_df.columns


@pytest.mark.parametrize("fname,sep,seq_col", [(f, s, c) for f, (s, c) in FILES.items()])
def test_mangled_casing(fname, sep, seq_col):
    """Randomly-cased headers (as if a user hand-edited or a tool changed
    version) must still map, now that matching is case-insensitive."""
    df = _read(fname, sep)

    def mangle(c):
        return "".join(ch.upper() if i % 2 == 0 else ch.lower() for i, ch in enumerate(c))

    df = df.rename(columns={c: mangle(c) for c in df.columns})
    content = _to_bytes(df, sep)
    result_df, status, err, warn = data_loader.load_and_validate_file(content, fname, "Peptidomic")
    assert status == "yes", f"{fname} (mangled headers): {err}"
    assert "Sequence" in result_df.columns
    assert "Protein" in result_df.columns
