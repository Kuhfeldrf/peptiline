"""
UniProt fetch and protein combination logic.
Extracted from notebook ProteinCombinationHandler class (Cell 3).
"""
import time
import re

import pandas as pd
import numpy as np

from .data_loader import find_species, extract_protein_id

from utils.uniprot_client import fetch_uniprot_info_batch, fetch_uniprot_info, UniProtClient


def extract_clean_protein_id(protein_string):
    """Extract clean protein ID from accession string."""
    if not protein_string or protein_string == 'Unknown':
        return protein_string
    protein_str = str(protein_string).strip()
    if '|' in protein_str:
        parts = protein_str.split('|')
        if len(parts) >= 3:
            # UniProt-style db|ACCESSION|NAME
            return parts[1].split(';')[0].split('/')[0].strip()
        elif len(parts) == 2:
            # PEAKS-style ACCESSION|NAME
            return parts[0].split(';')[0].split('/')[0].strip()
    return protein_str


def extract_protein_info_from_accession(protein_accession, protein_dict):
    """Extract protein name and species from UniProt accession format."""
    name = ""
    species = ""
    if "|" in protein_accession:
        parts = protein_accession.split("|")
        if len(parts) == 3:
            name_species = parts[2]
            if "_" in name_species:
                name, species = name_species.split("_", 1)
                species = find_species(species)
            else:
                name = name_species

    if name or species:
        protein_dict[protein_accession] = {
            "name": name if name else protein_accession,
            "species": species if species else "Unknown"
        }


def fetch_missing_proteins(missing_protein_ids, protein_dict, progress_callback=None):
    """
    Fetch missing proteins from UniProt.

    Args:
        missing_protein_ids: set of protein IDs to fetch
        protein_dict: dict to update with fetched info
        progress_callback: callable(current, total, message) for progress updates

    Returns:
        tuple: (success_count, dead_proteins set)
    """
    client = UniProtClient()
    protein_ids = list(missing_protein_ids)
    dead_proteins = set()
    bad_protein = set()

    batch_size = 10
    success_count = 0
    processed_count = 0
    total = len(protein_ids)

    # Process in batches
    batches = [protein_ids[i:i+batch_size] for i in range(0, total, batch_size)]

    for batch_idx, batch in enumerate(batches):
        if progress_callback:
            progress_callback(processed_count, total, f"Processing batch {batch_idx+1}/{len(batches)}")

        try:
            results = client.fetch_proteins_batch(batch)
            for protein_id in batch:
                if protein_id in results:
                    name, uniprot_species, signal_end = results[protein_id]
                    protein_dict[protein_id] = {
                        "name": name if name else protein_id,
                        "species": uniprot_species,
                        "signal_end": signal_end,
                    }
                    success_count += 1
                else:
                    bad_protein.add(protein_id)
                processed_count += 1
        except Exception:
            bad_protein.update(batch)
            processed_count += len(batch)

        time.sleep(0.5)

    # Individual fetch for failed proteins
    for protein_id in list(bad_protein):
        try:
            name, species, _, signal_end = fetch_uniprot_info(protein_id)
            if name is not None and species is not None:
                protein_dict[protein_id] = {"name": name, "species": species, "signal_end": signal_end}
                bad_protein.discard(protein_id)
                success_count += 1
            else:
                dead_proteins.add(protein_id)
                bad_protein.discard(protein_id)
                protein_dict[protein_id] = {"name": protein_id, "species": "Unknown", "signal_end": None}
        except Exception:
            dead_proteins.add(protein_id)
            bad_protein.discard(protein_id)
            protein_dict[protein_id] = {"name": protein_id, "species": "Unknown", "signal_end": None}
        time.sleep(0.5)

    if progress_callback:
        progress_callback(total, total, f"Complete: {success_count}/{total} proteins fetched")

    return success_count, dead_proteins


def extract_row_proteins(row):
    """Extract protein set from a row."""
    row_proteins = set()

    pos_value = row.get('Positions in Proteins', None)
    try:
        if pos_value is not None:
            pos_str = str(pos_value) if not pd.isna(pos_value) else None
            if pos_str and pos_str != 'Unknown' and pos_str != 'nan' and pos_str != 'None':
                if '; ' in pos_str:
                    for part in pos_str.split('; '):
                        if part.strip():
                            protein_id = part.strip().split()[0].split('[')[0]
                            if protein_id:
                                row_proteins.add(protein_id)
                elif '/ ' in pos_str:
                    for part in pos_str.split('/ '):
                        if part.strip():
                            protein_id = part.strip().split()[0].split('[')[0]
                            if protein_id:
                                row_proteins.add(protein_id)
                else:
                    parts = pos_str.strip().split()
                    if parts:
                        protein_id = parts[0].split('[')[0]
                        if protein_id:
                            row_proteins.add(protein_id)
    except (ValueError, TypeError):
        pass

    if not row_proteins:
        prot_value = row.get('Protein', row.get('Master Protein Accessions', None))
        try:
            if prot_value is None:
                return row_proteins
            try:
                if pd.isna(prot_value):
                    return row_proteins
            except (ValueError, TypeError):
                pass

            prot_str = str(prot_value)
            if prot_str and prot_str != 'Unknown' and prot_str != 'nan' and prot_str != 'None':
                if isinstance(prot_value, list):
                    if prot_value and any(p is not None and str(p).strip() and str(p) != 'Unknown' for p in prot_value):
                        row_proteins = set(extract_clean_protein_id(str(p).strip()) for p in prot_value if p)
                else:
                    if '; ' in prot_str:
                        row_proteins = set(extract_clean_protein_id(p.strip()) for p in prot_str.split('; ') if p.strip())
                    elif '/ ' in prot_str:
                        row_proteins = set(extract_clean_protein_id(p.strip()) for p in prot_str.split('/ ') if p.strip())
                    elif ',' in prot_str:
                        row_proteins = set(extract_clean_protein_id(p.strip()) for p in prot_str.split(',') if p.strip())
                    else:
                        clean_id = extract_clean_protein_id(prot_str.strip())
                        if clean_id:
                            row_proteins = {clean_id}
        except Exception:
            return set()

    return row_proteins


def get_protein_combinations(df, protein_dict):
    """Extract unique multi-protein combinations from the dataset."""
    if df is None or df.empty:
        return []

    protein_combinations = set()

    # Prepare working dataframe
    cols_needed = ["Positions in Proteins", "Protein"]
    if "Master Protein Accessions" in df.columns:
        working_df = df[["Positions in Proteins", "Master Protein Accessions"]].copy()
        working_df = working_df.rename(columns={"Master Protein Accessions": "Protein"})
    else:
        working_df = df.reindex(columns=cols_needed).copy()

    for col in cols_needed:
        if col not in working_df.columns:
            working_df[col] = np.nan

    working_df = working_df.fillna('Unknown')

    for row in working_df.itertuples(index=False, name=None):
        pos, prot = row
        try:
            position_proteins = []
            master_proteins = []

            if pos != 'Unknown':
                if '; ' in pos:
                    for part in pos.split('; '):
                        if part.strip():
                            protein_id = part.strip().split()[0].split('[')[0]
                            if protein_id:
                                position_proteins.append(protein_id)
                elif '/ ' in pos:
                    for part in pos.split('/ '):
                        if part.strip():
                            protein_id = part.strip().split()[0].split('[')[0]
                            if protein_id:
                                position_proteins.append(protein_id)
                else:
                    tokens = pos.strip().split()
                    if tokens:
                        protein_id = tokens[0].split('[')[0]
                        if protein_id:
                            position_proteins = [protein_id]

            if prot != 'Unknown':
                if isinstance(prot, list):
                    master_proteins = [extract_clean_protein_id(str(p).strip()) for p in prot if p]
                else:
                    master_str = str(prot)
                    if '; ' in master_str:
                        master_proteins = [extract_clean_protein_id(p.strip()) for p in master_str.split('; ') if p.strip()]
                    elif '/ ' in master_str:
                        master_proteins = [extract_clean_protein_id(p.strip()) for p in master_str.split('/ ') if p.strip()]
                    elif ',' in master_str:
                        master_proteins = [extract_clean_protein_id(p.strip()) for p in master_str.split(',') if p.strip()]
                    else:
                        clean_id = extract_clean_protein_id(master_str.strip())
                        if clean_id:
                            master_proteins = [clean_id]

            if len(position_proteins) > 1:
                protein_combinations.add('; '.join(sorted(position_proteins)))
            elif len(master_proteins) > 1:
                if position_proteins:
                    protein_combinations.add('; '.join(sorted(position_proteins)))
                else:
                    protein_combinations.add('; '.join(sorted(master_proteins)))
            elif pos == 'Unknown' and prot == 'Unknown':
                protein_combinations.add('Unknown')
            elif pos == 'Unknown' and len(master_proteins) > 1:
                protein_combinations.add('; '.join(sorted(master_proteins)))
            else:
                all_proteins = position_proteins + master_proteins
                species_set = set()
                for protein in all_proteins:
                    if protein in protein_dict:
                        # treat empty/falsy species same as Unknown
                        sp = protein_dict[protein].get('species') or 'Unknown'
                        species_set.add(sp)
                    elif protein == 'Unknown':
                        species_set.add('Unknown')
                    # proteins absent from protein_dict (but not the literal 'Unknown'
                    # string) are not added here — they appear in the missing-proteins
                    # warning instead of polluting the single-protein combinations list
                if len(species_set) > 1 or 'Unknown' in species_set:
                    if position_proteins:
                        protein_combinations.add('; '.join(sorted(position_proteins)))
                    elif master_proteins:
                        protein_combinations.add('; '.join(sorted(master_proteins)))

        except Exception:
            continue

    return list(protein_combinations)


def get_combination_details(combinations, df, protein_dict):
    """
    Get details for each protein combination for display.

    Returns list of dicts: [{
        'combo': 'P1; P2',
        'occurrences': 5,
        'proteins': [
            {'id': 'P1', 'species': 'cow', 'name': 'Alpha-casein', 'default_decision': 'new'},
            ...
        ]
    }, ...]
    """
    # Build lookup: combo_key -> list of row indices
    combo_to_rows = {}
    for idx in range(len(df)):
        try:
            row = df.iloc[idx]
            row_proteins = extract_row_proteins(row)
            if row_proteins:
                combo_key = '; '.join(sorted(row_proteins))
                if combo_key not in combo_to_rows:
                    combo_to_rows[combo_key] = []
                combo_to_rows[combo_key].append(idx)
        except Exception:
            continue

    # Pre-extract info for proteins not in dict
    all_proteins = set()
    for combo in combinations:
        if '; ' in combo:
            all_proteins.update(combo.split('; '))
        elif '/ ' in combo:
            all_proteins.update(combo.split('/ '))
        else:
            all_proteins.add(combo)

    for protein in all_proteins:
        if protein not in protein_dict and protein != 'Unknown':
            extract_protein_info_from_accession(protein, protein_dict)

    result = []
    for combo in combinations:
        if '; ' in combo:
            proteins = combo.split('; ')
        elif '/ ' in combo:
            proteins = combo.split('/ ')
        else:
            proteins = [combo]

        combo_row_indices = combo_to_rows.get(combo, [])
        combo_rows = [df.iloc[idx] for idx in combo_row_indices]
        occurrences = len(combo_rows)

        protein_details = []
        for protein in proteins:
            if protein == 'Unknown':
                species = "Unknown"
                name = "Unknown Protein"
            else:
                species = protein_dict.get(protein, {}).get('species', "Unknown")
                name = protein_dict.get(protein, {}).get('name', "Unknown Protein")

            default_decision = _get_default_decision(protein, combo_rows)

            protein_details.append({
                'id': protein,
                'species': species,
                'name': name,
                'default_decision': default_decision,
            })

        result.append({
            'combo': combo,
            'occurrences': occurrences,
            'proteins': protein_details,
        })

    return result


def _get_default_decision(protein, combo_rows):
    """Determine default decision for a protein."""
    if protein == 'Unknown':
        return 'asis'

    if not combo_rows:
        return 'new'

    first_row = combo_rows[0]

    if 'Leading razor protein' in first_row.index:
        leading_razor = first_row.get('Leading razor protein', '')
        if leading_razor:
            leading_razor_clean = extract_clean_protein_id(str(leading_razor))
            protein_clean = extract_clean_protein_id(protein)
            return 'new' if protein_clean == leading_razor_clean else 'remove'
        return 'new'

    master_proteins = []
    if isinstance(first_row.get('Protein'), list):
        if first_row['Protein'] and any(p is not None and str(p).strip() for p in first_row['Protein']):
            master_proteins = [str(p).strip() for p in first_row['Protein']]
    elif not pd.isna(first_row.get('Protein', pd.NA)):
        prot = str(first_row['Protein'])
        if ';' in prot:
            master_proteins = [p.strip() for p in prot.split(';')]
        elif '/' in prot:
            master_proteins = [p.strip() for p in prot.split('/')]
        elif ',' in prot:
            master_proteins = [p.strip() for p in prot.split(',')]
        else:
            master_proteins = [prot.strip()]

    if len(master_proteins) == 1:
        return 'new' if protein in master_proteins else 'remove'
    return 'new'


def apply_protein_decisions(df, decisions, protein_dict):
    """
    Apply user protein mapping decisions to the dataframe.

    Args:
        df: pandas DataFrame
        decisions: dict of {combo_string: {protein_id: decision_string}}
        protein_dict: protein dictionary

    Returns:
        tuple: (processed_df, error_message or None)
    """
    # Validation
    for combo, protein_decisions in decisions.items():
        has_asis = any(d.upper() == 'ASIS' for d in protein_decisions.values())
        has_other = any(d.upper() in ('NEW', 'REMOVE') or d.upper().startswith('CUSTOM:')
                        for d in protein_decisions.values())
        if has_asis and has_other:
            return df, f"Combination '{combo}': ASIS cannot be used with other decision types"

        for protein, decision in protein_decisions.items():
            d_upper = decision.upper()
            if d_upper not in ('NEW', 'REMOVE', 'ASIS') and not d_upper.startswith('CUSTOM:'):
                return df, f"Protein '{protein}': Invalid decision '{decision}'"
            if d_upper.startswith('CUSTOM:') and not decision.split(':', 1)[1].strip():
                return df, f"Protein '{protein}': CUSTOM requires a protein ID after the colon"

    processed_df = df.copy()

    # Build index: combo -> row indices
    combo_to_indices = {}
    for idx in range(len(processed_df)):
        row = processed_df.iloc[idx]
        row_proteins = extract_row_proteins(row)
        if row_proteins:
            combo_key = '; '.join(sorted(row_proteins))
            if combo_key not in combo_to_indices:
                combo_to_indices[combo_key] = []
            combo_to_indices[combo_key].append(idx)

    rows_to_remove = set()
    new_rows = []

    for combo, protein_decisions in decisions.items():
        matched_indices = combo_to_indices.get(combo, [])
        if not matched_indices:
            continue

        all_asis = all(d.upper() == 'ASIS' for d in protein_decisions.values())
        if all_asis:
            continue

        has_new = any(d.upper() == 'NEW' for d in protein_decisions.values())
        has_remove = any(d.upper() == 'REMOVE' for d in protein_decisions.values())
        has_custom = any(d.upper().startswith('CUSTOM:') for d in protein_decisions.values())

        for idx in matched_indices:
            row = processed_df.iloc[idx]

            if has_new:
                for protein, decision in protein_decisions.items():
                    if decision.upper() == 'NEW':
                        new_row = row.copy()
                        new_row['Protein'] = protein

                        original_positions = row.get('Positions in Proteins', 'Unknown')
                        if original_positions and str(original_positions) != 'Unknown' and not pd.isna(original_positions):
                            if '; ' in str(original_positions):
                                position_parts = str(original_positions).split('; ')
                            elif '/ ' in str(original_positions):
                                position_parts = str(original_positions).split('/ ')
                            else:
                                position_parts = [str(original_positions)]

                            position_found = False
                            for pos_part in position_parts:
                                if pos_part.strip().startswith(protein + ' ') or pos_part.strip().startswith(protein + '['):
                                    new_row['Positions in Proteins'] = pos_part.strip()
                                    position_found = True
                                    if '[' in pos_part and ']' in pos_part:
                                        try:
                                            range_part = pos_part[pos_part.find('[')+1:pos_part.find(']')]
                                            if '-' in range_part:
                                                start_pos, end_pos = range_part.split('-')
                                                new_row['start'] = int(start_pos.strip())
                                                new_row['end'] = int(end_pos.strip())
                                        except Exception:
                                            pass
                                    break

                            if not position_found:
                                new_row['Positions in Proteins'] = 'Unknown'
                                if 'start' in new_row:
                                    new_row['start'] = np.nan
                                if 'end' in new_row:
                                    new_row['end'] = np.nan
                        else:
                            new_row['Positions in Proteins'] = 'Unknown'
                            if 'start' in new_row:
                                new_row['start'] = np.nan
                            if 'end' in new_row:
                                new_row['end'] = np.nan

                        new_rows.append(new_row)

                rows_to_remove.add(idx)

            elif has_custom or (has_remove and has_custom):
                proteins_in_combo = set(combo.split('; ')) if '; ' in combo else {combo}
                modified_proteins = proteins_in_combo.copy()
                custom_mappings = {}  # Track old->new protein ID mappings

                for protein, decision in protein_decisions.items():
                    d_upper = decision.upper()
                    if d_upper.startswith('CUSTOM:'):
                        new_protein_id = decision.split(':', 1)[1].strip()
                        custom_mappings[protein] = new_protein_id
                        if protein in modified_proteins:
                            modified_proteins.remove(protein)
                            modified_proteins.add(new_protein_id)
                    elif d_upper == 'REMOVE':
                        modified_proteins.discard(protein)

                if modified_proteins:
                    processed_df.at[processed_df.index[idx], 'Protein'] = '; '.join(sorted(modified_proteins))

                    # Update Positions in Proteins with custom mappings, mirroring notebook logic
                    positions = row.get('Positions in Proteins', '')
                    if positions and str(positions) != 'Unknown' and not pd.isna(positions):
                        new_positions = []
                        updated_start = None
                        updated_end = None

                        pos_str = str(positions)
                        if '; ' in pos_str:
                            pos_parts = pos_str.split('; ')
                        elif '/ ' in pos_str:
                            pos_parts = pos_str.split('/ ')
                        else:
                            pos_parts = [pos_str]

                        for pos_part in pos_parts:
                            pos_modified = False

                            # Replace old protein ID with new one for CUSTOM mappings
                            for old_protein, new_protein in custom_mappings.items():
                                if pos_part.strip().startswith(old_protein + ' ') or pos_part.strip().startswith(old_protein + '['):
                                    if '[' in pos_part and ']' in pos_part:
                                        range_start_idx = pos_part.find('[')
                                        range_end_idx = pos_part.find(']') + 1
                                        range_part = pos_part[range_start_idx:range_end_idx]
                                        new_pos = f"{new_protein} {range_part}"
                                        # Extract start/end when result is a single protein
                                        try:
                                            range_content = range_part[1:-1]  # strip brackets
                                            if '-' in range_content:
                                                start_str, end_str = range_content.split('-', 1)
                                                if len(modified_proteins) == 1:
                                                    updated_start = int(start_str.strip())
                                                    updated_end = int(end_str.strip())
                                        except Exception:
                                            pass
                                    else:
                                        new_pos = new_protein
                                    new_positions.append(new_pos)
                                    pos_modified = True
                                    break

                            # If not a custom-mapped position, keep it if its protein is still present
                            if not pos_modified:
                                for protein in modified_proteins:
                                    if pos_part.strip().startswith(protein + ' ') or pos_part.strip().startswith(protein + '['):
                                        new_positions.append(pos_part.strip())
                                        # Extract start/end when result is a single protein
                                        if len(modified_proteins) == 1 and '[' in pos_part and ']' in pos_part:
                                            try:
                                                range_part = pos_part[pos_part.find('[') + 1:pos_part.find(']')]
                                                if '-' in range_part:
                                                    start_str, end_str = range_part.split('-', 1)
                                                    updated_start = int(start_str.strip())
                                                    updated_end = int(end_str.strip())
                                            except Exception:
                                                pass
                                        break

                        if new_positions:
                            processed_df.at[processed_df.index[idx], 'Positions in Proteins'] = '; '.join(new_positions)
                            # Update start/end columns for single-protein result
                            if updated_start is not None and updated_end is not None:
                                if 'start' in processed_df.columns:
                                    processed_df.at[processed_df.index[idx], 'start'] = updated_start
                                if 'end' in processed_df.columns:
                                    processed_df.at[processed_df.index[idx], 'end'] = updated_end
                        else:
                            processed_df.at[processed_df.index[idx], 'Positions in Proteins'] = 'Unknown'
                else:
                    rows_to_remove.add(idx)

            elif has_remove and not has_custom:
                proteins_in_combo = set(combo.split('; ')) if '; ' in combo else {combo}
                proteins_to_remove = {p for p, d in protein_decisions.items() if d.upper() == 'REMOVE'}
                remaining = proteins_in_combo - proteins_to_remove

                if remaining:
                    processed_df.at[processed_df.index[idx], 'Protein'] = '; '.join(sorted(remaining))
                else:
                    rows_to_remove.add(idx)

    # Apply changes
    if rows_to_remove:
        indices_to_remove = [processed_df.index[i] for i in sorted(rows_to_remove, reverse=True)]
        processed_df = processed_df.drop(index=indices_to_remove)

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        processed_df = pd.concat([processed_df, new_df], ignore_index=True)

    return processed_df, None


# ---------------------------------------------------------------------------
# Merge / rename protein sources
# ---------------------------------------------------------------------------
#
# The combination controls above only act on multi-protein rows.  Single-protein
# rows (e.g. a peptide that only ever maps to the beta-casein A2 variant
# ``P02666A2``) never surface there, so they cannot be folded into a canonical
# accession from the UI.  These helpers add an accession-level rename/merge that
# operates on *every* row and rewrites only the derived columns
# (``Protein``/``protein_name``/``protein_species``/``Positions in Proteins``),
# deliberately leaving the original provenance columns
# (``Master Protein Accessions`` and the bioactive-search ``protein_id`` /
# ``protein_description`` / ``species``) untouched.


def _split_protein_members(value):
    """Split a Protein cell into an ordered list of clean member accessions."""
    if value is None:
        return []
    try:
        if isinstance(value, list):
            return [extract_clean_protein_id(str(p).strip()) for p in value
                    if p is not None and str(p).strip()]
        if pd.isna(value):
            return []
    except (ValueError, TypeError):
        pass

    text = str(value).strip()
    if not text or text in ('Unknown', 'nan', 'None'):
        return []

    if '; ' in text:
        parts = text.split('; ')
    elif '/ ' in text:
        parts = text.split('/ ')
    elif ',' in text:
        parts = text.split(',')
    else:
        parts = [text]
    return [extract_clean_protein_id(p.strip()) for p in parts if p.strip()]


def _split_position_parts(value):
    """Split a 'Positions in Proteins' cell into its individual '<ID> [range]' parts."""
    if value is None:
        return []
    try:
        if pd.isna(value):
            return []
    except (ValueError, TypeError):
        pass
    text = str(value).strip()
    if not text or text in ('Unknown', 'nan', 'None'):
        return []
    if '; ' in text:
        return [p.strip() for p in text.split('; ') if p.strip()]
    if '/ ' in text:
        return [p.strip() for p in text.split('/ ') if p.strip()]
    return [text]


def get_protein_sources(df, protein_dict):
    """
    List every distinct protein accession present in the data for the
    merge/rename picker.

    Returns a list of dicts sorted by id:
        [{'id', 'name', 'species', 'count'}, ...]
    where ``count`` is the number of rows the accession appears in (whether as a
    single-protein row or as one member of a multi-protein row).
    """
    if df is None or getattr(df, 'empty', True):
        return []

    counts = {}
    for idx in range(len(df)):
        try:
            row_proteins = extract_row_proteins(df.iloc[idx])
        except Exception:
            continue
        for pid in row_proteins:
            if pid and pid != 'Unknown':
                counts[pid] = counts.get(pid, 0) + 1

    sources = []
    for pid in sorted(counts):
        info = protein_dict.get(pid, {})
        sources.append({
            'id': pid,
            'name': info.get('name', ''),
            'species': info.get('species', ''),
            'count': counts[pid],
        })
    return sources


def _normalize_rename_groups(groups):
    """Validate and normalize merge/rename groups. Returns (clean_groups, error)."""
    clean = []
    seen_sources = {}
    for i, group in enumerate(groups or []):
        if not isinstance(group, dict):
            return None, f"Merge group {i + 1}: invalid format"
        sources = [str(s).strip() for s in group.get('sources', []) if str(s).strip()]
        target_id = str(group.get('target_id', '')).strip()
        if not sources:
            return None, f"Merge group {i + 1}: select at least one protein source"
        if not target_id:
            return None, f"Merge group {i + 1}: a target protein ID is required"
        for s in sources:
            if s in seen_sources and seen_sources[s] != target_id:
                return None, (f"Protein '{s}' is assigned to two different targets "
                              f"('{seen_sources[s]}' and '{target_id}')")
            seen_sources[s] = target_id
        clean.append({
            'sources': sources,
            'target_id': target_id,
            'target_name': str(group.get('target_name', '')).strip(),
        })
    return clean, None


def apply_source_renames(df, groups, protein_dict):
    """
    Merge/rename protein accessions across every row.

    Args:
        df: pandas DataFrame
        groups: list of {'sources': [id...], 'target_id', 'target_name'}
        protein_dict: protein dictionary (updated in place for target accessions)

    Returns:
        tuple: (processed_df, error_message or None)
    """
    clean_groups, error = _normalize_rename_groups(groups)
    if error:
        return df, error
    if not clean_groups:
        return df, None

    # source id -> target id, and target id -> resolved name
    source_to_target = {}
    target_info = {}
    for group in clean_groups:
        tid = group['target_id']
        for s in group['sources']:
            source_to_target[s] = tid
        # Later groups win if a target is defined twice; only overwrite with
        # a non-empty override so a blank field doesn't clobber a real value.
        info = target_info.setdefault(tid, {'name': ''})
        if group['target_name']:
            info['name'] = group['target_name']

    # Fill any unspecified name from the existing protein_dict entry. Species
    # is intentionally not tracked by merge/rename groups; existing per-row
    # species values (and protein_dict species) are left untouched.
    for tid, info in target_info.items():
        if not info['name']:
            info['name'] = protein_dict.get(tid, {}).get('name', '') or ''
        # Keep protein_dict consistent so add_protein_info fills 'Unknown' rows.
        if info['name']:
            protein_dict.setdefault(tid, {})['name'] = info['name']

    processed_df = df.copy()
    has_protein = 'Protein' in processed_df.columns
    has_positions = 'Positions in Proteins' in processed_df.columns
    has_name = 'protein_name' in processed_df.columns

    for idx in range(len(processed_df)):
        row = processed_df.iloc[idx]
        row_label = processed_df.index[idx]

        # --- Protein column ---
        new_single_target = None
        if has_protein:
            members = _split_protein_members(row.get('Protein'))
            if members:
                mapped = [source_to_target.get(m, m) for m in members]
                # de-duplicate, preserve nothing fancy — join sorted like the
                # rest of this module does.
                unique_mapped = sorted(set(mapped))
                if unique_mapped != sorted(set(members)):
                    processed_df.at[row_label, 'Protein'] = '; '.join(unique_mapped)
                if len(unique_mapped) == 1 and unique_mapped[0] in target_info:
                    new_single_target = unique_mapped[0]

        # --- Positions in Proteins column ---
        if has_positions:
            parts = _split_position_parts(row.get('Positions in Proteins'))
            if parts:
                new_parts = []
                changed = False
                for part in parts:
                    token = part.split()[0].split('[')[0] if part.split() else ''
                    if token in source_to_target:
                        tid = source_to_target[token]
                        remainder = part[len(token):]
                        new_part = tid + remainder
                        changed = True
                    else:
                        new_part = part
                    if new_part not in new_parts:
                        new_parts.append(new_part)
                if changed:
                    processed_df.at[row_label, 'Positions in Proteins'] = '; '.join(new_parts)

        # --- Unify name for rows that resolve to a single target ---
        if new_single_target is not None:
            info = target_info[new_single_target]
            if has_name and info['name']:
                processed_df.at[row_label, 'protein_name'] = info['name']

    return processed_df, None


def detect_rename_conflicts(source_renames, protein_decisions):
    """
    Warn when a merge/rename group touches an accession that a combination
    decision already acted on. The merge/rename is authoritative and still
    applies; these are non-fatal advisories so the user knows the earlier
    combined-protein choice is being overridden.

    Args:
        source_renames: list of merge groups (as sent by the UI)
        protein_decisions: saved combination decisions, either simplified UI
            format ({combo: {action, ...}}) or backend format
            ({combo: {protein_id: decision_str}})

    Returns:
        list of human-readable warning strings.
    """
    clean_groups, error = _normalize_rename_groups(source_renames)
    if error or not clean_groups:
        return []

    # combo string -> set of accessions the combination step acted on (non-ASIS)
    touched = {}
    for combo, data in (protein_decisions or {}).items():
        if not isinstance(data, dict):
            continue
        combo_ids = set()
        if 'action' in data:
            action = str(data.get('action', 'ASIS')).upper()
            if action == 'ASIS':
                continue
            combo_ids.update(p.strip() for p in combo.split(';') if p.strip())
            pid = str(data.get('protein_id', '')).strip()
            if pid:
                combo_ids.add(pid)
            combo_ids.update(str(p).strip() for p in data.get('protein_ids', []))
        else:
            # backend format {protein_id: decision_str}
            non_asis = False
            for pid, dec in data.items():
                dec_str = str(dec).upper()
                if dec_str == 'ASIS':
                    continue
                non_asis = True
                combo_ids.add(str(pid).strip())
                if dec_str.startswith('CUSTOM:'):
                    combo_ids.add(dec_str.split('CUSTOM:', 1)[1].strip())
            if not non_asis:
                continue
        if combo_ids:
            touched[combo] = combo_ids

    warnings = []
    for group in clean_groups:
        group_ids = set(group['sources']) | {group['target_id']}
        for combo, combo_ids in touched.items():
            overlap = sorted(group_ids & combo_ids)
            if overlap:
                warnings.append(
                    f"Merge/rename of {', '.join(overlap)} overrides your "
                    f"combined-protein decision for '{combo}'."
                )
    return warnings


# ---------------------------------------------------------------------------
# Protein-mapping-key JSON format: fold merge/rename groups in and out of the
# same ``protein_decisions`` dict the combination step already saves, using
# ``"action": "MERGE"`` entries, so the downloadable/uploadable mapping key
# has a single unified shape instead of a second top-level array.
# ---------------------------------------------------------------------------

def merge_groups_as_decisions(groups):
    """
    Represent merge/rename groups as ``protein_decisions``-style entries.

    Each group becomes one entry keyed by its sorted, semicolon-joined source
    accessions (matching the combo-key convention already used for
    ASIS/SPLIT/CUSTOM entries), with ``action: "MERGE"``.
    """
    clean_groups, error = _normalize_rename_groups(groups)
    if error or not clean_groups:
        return {}
    entries = {}
    for group in clean_groups:
        key = '; '.join(sorted(group['sources']))
        entries[key] = {
            'action': 'MERGE',
            'protein_ids': group['sources'],
            'target_id': group['target_id'],
            'target_name': group['target_name'],
        }
    return entries


def extract_merge_groups(protein_decisions):
    """
    Inverse of ``merge_groups_as_decisions``: pull any ``action: "MERGE"``
    entries out of an uploaded ``protein_decisions`` dict.

    Returns (remaining_decisions, groups) where ``remaining_decisions`` has
    the MERGE entries removed (safe to hand to the combination-decision
    pipeline) and ``groups`` is the internal merge/rename group format
    (``{'sources', 'target_id', 'target_name'}``).
    """
    remaining = {}
    groups = []
    for combo, data in (protein_decisions or {}).items():
        if isinstance(data, dict) and str(data.get('action', '')).upper() == 'MERGE':
            sources = [str(s).strip() for s in data.get('protein_ids') or combo.split('; ')
                       if str(s).strip()]
            groups.append({
                'sources': sources,
                'target_id': str(data.get('target_id', '')).strip(),
                'target_name': str(data.get('target_name', '')).strip(),
            })
        else:
            remaining[combo] = data
    return remaining, groups
