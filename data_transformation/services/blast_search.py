"""
BLAST search logic and database queries.
Extracted from notebook DataTransformation class (Cell 1).
Converted from raw sqlite3 to Django ORM.
"""
import os
import csv
import subprocess
import time
import shutil
from collections import defaultdict

import pandas as pd
import numpy as np
from django.conf import settings

# Queries run against mbpdb_replica's models -- a read-only copy of MBPDB's
# ProteinInfo/PeptideInfo/Function/Reference tables. mbpdb_replica.sqlite3 is
# committed to the repo already populated; `manage.py loadreplica <db.sqlite3>`
# refreshes it from a newer MBPDB dump. If the replica is ever empty a search
# returns no results rather than crashing the Data Transformation dashboard.
from mbpdb_replica.models import PeptideInfo, ProteinInfo, Function, Reference


def create_work_directory(base_dir=None):
    """Create a working directory for BLAST operations."""
    if base_dir is None:
        base_dir = settings.WORK_DIRECTORY
    path = os.path.join(base_dir, f'work_{int(round(time.time() * 1000))}')
    os.makedirs(path, exist_ok=True)
    return path


def cleanup_work_directories(work_directory=None, keep=12, max_age_hours=8):
    """Clean up old wizard work directories.

    Two rules, applied in order:

    * age — anything older than ``max_age_hours`` is removed unconditionally.
      Sessions time out at 4 h, so 8 h is a orphaned-for-certain cutoff and
      keeps the temp dir from growing without bound between deploys.
    * count — of what remains, only the ``keep`` most-recent survive.

    The per-search BLAST scratch (replica FASTA, makeblastdb index, query file,
    tabular output) is *not* covered here: run_blast_search removes its own
    ``_blast_scratch`` directory in a ``finally`` the moment the search ends, so
    those never accumulate inside a live session dir in the first place.
    """
    if work_directory is None:
        work_directory = settings.WORK_DIRECTORY
    try:
        dirs = [f for f in os.scandir(work_directory) if f.is_dir() and f.name.startswith('work_')]
        cutoff = time.time() - max_age_hours * 3600
        survivors = []
        for entry in dirs:
            try:
                if os.path.getmtime(entry.path) < cutoff:
                    shutil.rmtree(entry.path, ignore_errors=True)
                else:
                    survivors.append(entry)
            except Exception:
                pass

        survivors.sort(key=lambda x: os.path.getmtime(x.path), reverse=True)
        for dir_entry in survivors[keep:]:
            shutil.rmtree(dir_entry.path, ignore_errors=True)
    except Exception:
        pass


def make_blast_db(library_fasta_path):
    """Create BLAST database from FASTA file."""
    subprocess.check_output(
        ['makeblastdb', '-in', library_fasta_path, '-dbtype', 'prot'],
        stderr=subprocess.STDOUT
    )


def run_blast_search(peptides, similarity_threshold=100, work_dir=None,
                     progress_callback=None):
    """
    Search peptides against the MBPDB replica.

    Two passes:

    * Exact match (threshold 100, or peptides shorter than the BLAST word size)
      -- a SINGLE ``peptide__in=[...]`` query with ``functions__references``
      prefetched, instead of one query per peptide plus a Function/Reference
      N+1. No BLAST database is built for this pass.
    * BLAST (threshold < 100) -- all remaining peptides go through one batched
      ``blastp`` call; every subject hit across every query is then resolved in a
      SINGLE ``id__in=[...]`` query. All BLAST scratch (the replica FASTA, the
      makeblastdb index, the query file, the tabular output) is written under a
      per-search ``_blast_scratch`` directory and removed in a ``finally`` so it
      never accumulates.

    Args:
        peptides: list of peptide sequences
        similarity_threshold: minimum similarity percentage (0-100)
        work_dir: working directory for temp files
        progress_callback: callable(current, total, status_msg) for progress updates

    Returns:
        pd.DataFrame with search results
    """
    if work_dir is None:
        work_dir = create_work_directory()

    total = len(peptides)
    results = []

    # Partition into exact-match vs BLAST peptides.
    exact_seqs = set()
    blast_peptides = []          # [(original_index, peptide)]
    for idx, peptide in enumerate(peptides):
        if similarity_threshold == 100 or len(peptide) < 4:
            exact_seqs.add(peptide)
        else:
            blast_peptides.append((idx, peptide))

    # --- Exact-match pass: ONE query for every distinct exact sequence ---
    if exact_seqs:
        if progress_callback:
            progress_callback(0, total, f"Exact match ({len(exact_seqs)} sequences)")
        exact_qs = (PeptideInfo.objects
                    .filter(peptide__in=exact_seqs)
                    .select_related('protein')
                    .prefetch_related('functions__references'))
        exact_rows = []
        for p in exact_qs:
            exact_rows.extend(_peptide_rows(p, p.peptide))
        if exact_rows:
            results.append(pd.DataFrame(exact_rows))

    # --- BLAST pass: one subprocess for all peptides; scratch removed after ---
    if blast_peptides:
        blast_dir = os.path.join(work_dir, '_blast_scratch')
        os.makedirs(blast_dir, exist_ok=True)
        try:
            if progress_callback:
                progress_callback(len(exact_seqs), total,
                                  f"Running BLAST on {len(blast_peptides)} peptides...")

            # Replica FASTA + makeblastdb index -- only built when BLAST is
            # actually needed (a pure exact-match search skips this entirely).
            fasta_db_path = os.path.join(blast_dir, "db.fasta")
            db_rows = 0
            with open(fasta_db_path, 'w') as f:
                for pep_id, pep_seq in PeptideInfo.objects.values_list('id', 'peptide').iterator():
                    f.write(f">{pep_id}\n{pep_seq}\n")
                    db_rows += 1

            if db_rows == 0:
                # An empty replica -- almost always mbpdb_replica.sqlite3 not
                # shipped in the deploy image (it's excluded from .dockerignore's
                # copy, or a stale build). makeblastdb on an empty FASTA would
                # otherwise raise a cryptic CalledProcessError; surface the real
                # cause instead.
                raise RuntimeError(
                    "The MBPDB reference database (mbpdb_replica) is empty -- no "
                    "peptides to search against. The bundled snapshot is likely "
                    "missing from this deployment. Upload a functional annotation "
                    "file instead, or contact the site administrator."
                )

            make_blast_db(fasta_db_path)

            # Multi-query FASTA; query IDs "q0", "q1", ... map hits back to peptides.
            query_path = os.path.join(blast_dir, "query.fasta")
            query_id_to_peptide = {}
            with open(query_path, 'w') as qf:
                for i, (_orig_idx, peptide) in enumerate(blast_peptides):
                    query_id = f"q{i}"
                    query_id_to_peptide[query_id] = peptide
                    qf.write(f">{query_id}\n{peptide}\n")

            output_path = os.path.join(blast_dir, "blastp_short.out")
            blast_args = [
                "blastp",
                "-query", query_path,
                "-db", fasta_db_path,
                "-outfmt", "6 std ppos qcovs qlen slen positive",
                "-evalue", "1000",
                "-word_size", "2",
                "-matrix", "IDENTITY",
                "-threshold", "1",
                "-task", "blastp-short",
                "-out", output_path,
            ]
            try:
                subprocess.check_output(blast_args, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                query_hits = {}     # no hits / BLAST error -- keep exact results
            else:
                query_hits = _process_blast_results_multi(output_path, similarity_threshold)

            if query_hits:
                # ONE query for every subject hit across ALL query peptides
                # (was one query per hitting peptide, each with its own N+1).
                all_subject_ids = set()
                for hits in query_hits.values():
                    all_subject_ids.update(hits)
                int_ids = [i for i in (_as_int(s) for s in all_subject_ids) if i is not None]
                subj_by_id = {
                    p.id: p
                    for p in (PeptideInfo.objects
                              .filter(id__in=int_ids)
                              .select_related('protein')
                              .prefetch_related('functions__references'))
                }
                for query_id, peptide in query_id_to_peptide.items():
                    rows = []
                    for sid, detail in query_hits.get(query_id, {}).items():
                        p = subj_by_id.get(_as_int(sid))
                        if p is not None:
                            rows.extend(_peptide_rows(p, peptide, blast_detail=detail))
                    if rows:
                        results.append(pd.DataFrame(rows))
        finally:
            shutil.rmtree(blast_dir, ignore_errors=True)

    if progress_callback:
        progress_callback(total, total, "Search complete")

    cleanup_work_directories()

    return _combine_results(results)


def _as_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


# BLAST tabular detail columns, in the order _process_blast_results_multi packs them.
_BLAST_DETAIL_COLS = (
    '% Alignment', 'Query start', 'Query end', 'Subject start', 'Subject end',
    'e-value', 'Alignment length', 'Mismatches', 'Gap opens',
)

# The MBPDB-shaped columns every result row carries.
_MBPDB_COLUMNS = [
    'search_peptide', 'protein_id', 'peptide', 'protein_description',
    'species', 'intervals', 'function', 'additional_details', 'ic50',
    'inhibition_type', 'inhibited_microorganisms', 'ptm', 'title',
    'authors', 'abstract', 'doi', 'search_type', 'scoring_matrix',
]


def _peptide_rows(pep_info, search_peptide, blast_detail=None):
    """Result rows for one ``PeptideInfo`` -- one per (function, reference), or a
    single ``function=None`` row when the peptide has no functions.

    ``pep_info`` must come from a queryset with ``select_related('protein')`` and
    ``prefetch_related('functions__references')`` so this makes no DB queries.
    ``blast_detail`` (optional) is the list from _process_blast_results_multi,
    merged in as the BLAST tabular columns.
    """
    base = {
        'search_peptide': search_peptide,
        'protein_id': pep_info.protein.pid,
        'peptide_id': pep_info.id,
        'peptide': pep_info.peptide,
        'protein_description': pep_info.protein.desc,
        'species': pep_info.protein.species,
        'intervals': pep_info.intervals,
        'search_type': 'sequence',
        'scoring_matrix': 'IDENTITY',
    }
    if blast_detail is not None:
        base.update(dict(zip(_BLAST_DETAIL_COLS, blast_detail)))

    funcs = list(pep_info.functions.all())
    if not funcs:
        return [{
            **base, 'function': None, 'additional_details': None, 'ic50': None,
            'inhibition_type': None, 'inhibited_microorganisms': None, 'ptm': None,
            'title': None, 'authors': None, 'abstract': None, 'doi': None,
        }]

    rows = []
    for func in funcs:
        for ref in func.references.all():
            rows.append({
                **base,
                'function': func.function,
                'additional_details': ref.additional_details,
                'ic50': ref.ic50,
                'inhibition_type': ref.inhibition_type,
                'inhibited_microorganisms': ref.inhibited_microorganisms,
                'ptm': ref.ptm,
                'title': ref.title,
                'authors': ref.authors,
                'abstract': ref.abstract,
                'doi': ref.doi,
            })
    return rows


def _process_blast_results_multi(output_path, similarity_threshold):
    """
    Parse a multi-query BLAST tabular output file.

    Returns a dict: query_id -> {subject_id: extra_info_list}
    """
    query_hits = defaultdict(dict)
    csv.register_dialect('blast_dialect', delimiter='\t')

    try:
        with open(output_path, 'r') as output_file:
            blast_data = csv.DictReader(
                output_file,
                fieldnames=['query', 'subject', 'percid', 'align_len', 'mismatches',
                            'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue',
                            'bitscore', 'ppos', 'qcov', 'qlen', 'slen', 'numpos'],
                dialect='blast_dialect'
            )

            for row in blast_data:
                tlen = max(float(row['slen']), float(row['qlen']))
                simcalc = 100 * ((float(row['numpos']) - float(row['gaps'])) / tlen)

                if simcalc >= similarity_threshold:
                    query_hits[row['query']][row['subject']] = [
                        f"{simcalc:.2f}", row['qstart'], row['qend'], row['sstart'],
                        row['send'], row['evalue'], row['align_len'], row['mismatches'],
                        row['gaps']
                    ]
    except FileNotFoundError:
        pass

    return query_hits


def _combine_results(results):
    """Combine and format final results."""
    if not results:
        return pd.DataFrame(columns=_MBPDB_COLUMNS)

    final_results = pd.concat(results, ignore_index=True)

    if 'peptide_id' in final_results.columns:
        final_results = final_results.drop('peptide_id', axis=1)

    sort_columns = ['search_peptide']
    if '% Alignment' in final_results.columns:
        sort_columns.append('% Alignment')

    return final_results.sort_values(
        sort_columns,
        ascending=[True] + [False] * (len(sort_columns) - 1)
    )


def format_search_results(final_results):
    """Format search results by aggregating duplicate groups."""
    if final_results.empty:
        return final_results

    if '% Alignment' in final_results.columns:
        final_results['% Alignment'] = pd.to_numeric(
            final_results['% Alignment'], errors='coerce'
        )

    grouped = final_results.groupby(["search_peptide", "function"], as_index=False)
    aggregated_results = []
    processed_indices = set()

    for _, group in grouped:
        if len(group) > 1:
            aggregated_row = _aggregate_group_data(group)
            aggregated_results.append(aggregated_row)
            processed_indices.update(group.index)

    remaining_rows = final_results.loc[~final_results.index.isin(processed_indices)]
    aggregated_df = pd.DataFrame(aggregated_results)

    return pd.concat([aggregated_df, remaining_rows], ignore_index=True)


def _aggregate_group_data(group):
    """Aggregate data for a group of results."""
    def enumerate_field(field):
        if field in group.columns and not group[field].dropna().empty:
            valid_values = set(group[field].dropna().astype(str).str.strip())
            valid_values = {val for val in valid_values if val != ''}
            if len(valid_values) > 1:
                return "; ".join([f"{i + 1}) {val}" for i, val in enumerate(valid_values)])
            elif len(valid_values) == 1:
                return next(iter(valid_values))
            return ''
        return ''

    return {col: enumerate_field(col) for col in group.columns}
