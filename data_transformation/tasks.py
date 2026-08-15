"""
Celery tasks for long-running data transformation operations.
"""
import re
import time
import os

from celery import shared_task
from django.core.cache import cache
from django.conf import settings

from .services.blast_search import run_blast_search, create_work_directory


@shared_task(bind=True)
def run_blast_search_task(self, work_dir, peptides, threshold):
    """
    Run BLAST search as a Celery task with progress tracking.

    Args:
        work_dir: working directory for temp files
        peptides: list of peptide sequences
        threshold: similarity threshold (0-100)

    Returns:
        dict with 'work_dir' and 'result_path' (parquet file path)
    """
    import pandas as pd

    task_id = self.request.id
    total = len(peptides)
    cache.set(f'size_{task_id}', total)
    start_time = time.time()

    def progress_callback(current, total_count, message):
        elapsed = time.time() - start_time
        cache.set(f'progress_{task_id}', current)
        cache.set(f'elapsed_time_{task_id}', elapsed)

    try:
        results_df = run_blast_search(
            peptides,
            similarity_threshold=threshold,
            work_dir=work_dir,
            progress_callback=progress_callback
        )

        result_path = os.path.join(work_dir, 'mbpdb_results.pkl')
        results_df.to_pickle(result_path)

        cache.set(f'progress_{task_id}', total)
        cache.set(f'elapsed_time_{task_id}', time.time() - start_time)
        cache.set(f'status_{task_id}', 'complete')

        return {
            'work_dir': work_dir,
            'result_path': result_path,
            'count': len(results_df),
        }
    except Exception as e:
        cache.set(f'status_{task_id}', 'failed')
        cache.set(f'error_{task_id}', str(e))
        raise


@shared_task(bind=True)
def fetch_uniprot_task(self, missing_protein_ids):
    """
    Fetch protein info from UniProt in batch as a Celery task.

    Args:
        missing_protein_ids: list of protein ID strings

    Returns:
        dict mapping protein_id -> {name, species, sequence, ...}
    """
    task_id = self.request.id
    total = len(missing_protein_ids)
    cache.set(f'size_{task_id}', total)
    start_time = time.time()

    from utils.uniprot_client import UniProtClient

    # UniProt accession IDs are 6 or 10 alphanumeric chars.
    # IDs like "P02666A1" (isoform suffixes, 7-9 chars) are invalid and cause
    # HTTP 400 for the entire batch, so filter them out before querying.
    # Valid 6-char:  [A-Z][0-9][A-Z0-9]{3}[0-9]  (e.g. P02666, Q9UKV3)
    # Valid 10-char: [A-Z][0-9][A-Z0-9]{7}[0-9]  (e.g. A0A023GPI8)
    _UNIPROT_ACC_RE = re.compile(
        r'^[A-Z][0-9][A-Z0-9]{3}[0-9]([A-Z0-9]{4})?$'
    )

    def _is_valid_uniprot_id(pid):
        return bool(_UNIPROT_ACC_RE.match(pid.upper())) if pid else False

    client = UniProtClient()
    found = {}
    BATCH_SIZE = 20

    valid_ids = [pid for pid in missing_protein_ids if _is_valid_uniprot_id(pid)]
    skipped_ids = [pid for pid in missing_protein_ids if not _is_valid_uniprot_id(pid)]
    if skipped_ids:
        print(f'fetch_uniprot_task: skipping {len(skipped_ids)} non-standard IDs: {skipped_ids[:10]}')

    # Use valid_ids count as the progress total so the bar reaches 100 %.
    cache.set(f'size_{task_id}', max(len(valid_ids), 1))

    processed = 0
    for batch_start in range(0, len(valid_ids), BATCH_SIZE):
        batch = valid_ids[batch_start:batch_start + BATCH_SIZE]

        elapsed = time.time() - start_time
        cache.set(f'progress_{task_id}', processed)
        cache.set(f'elapsed_time_{task_id}', elapsed)

        try:
            batch_results = client.fetch_proteins_batch(batch)
            for pid, info_tuple in batch_results.items():
                if info_tuple and len(info_tuple) >= 2:
                    name, species = info_tuple[0], info_tuple[1]
                    signal_end = info_tuple[2] if len(info_tuple) >= 3 else None
                    found[pid] = {'name': name, 'species': species, 'signal_end': signal_end}
        except Exception as exc:
            print(f'fetch_uniprot_task: batch error: {exc}')

        processed += len(batch)
        cache.set(f'progress_{task_id}', processed)

    # IDs that were valid UniProt format but returned no result
    not_found_ids = [pid for pid in valid_ids if pid not in found]

    cache.set(f'progress_{task_id}', max(len(valid_ids), 1))
    cache.set(f'elapsed_time_{task_id}', time.time() - start_time)
    cache.set(f'status_{task_id}', 'complete')

    # Store found proteins in Django cache as a fallback for views that cannot
    # retrieve the Celery result from the result backend (e.g. Redis eviction or
    # a race between the worker finishing and the web process polling).
    if found:
        cache.set(f'uniprot_found_{task_id}', found, 3600)

    return {
        'found': found,
        'not_found': not_found_ids,
        'skipped': skipped_ids,
    }
