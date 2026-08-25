"""
Refresh the mbpdb_replica DB alias by copying the five reference tables
straight out of a copy of MBPDB's own db.sqlite3, via SQLite's ATTACH
DATABASE -- no Django serialization round-trip, so the replica is a
byte-identical copy of MBPDB's data (same IDs, same values) and PeptiLine's
BLAST search returns exactly what MBPDB's own search would for the same
query. See docs/SPLIT_PLAN.md section 3.

MBPDB is the sole writer; PeptiLine never modifies these tables. Every run
wipes and reloads all five, so there's no merge/conflict logic to worry
about -- just "replace with the latest snapshot from the source file."
"""
from django.core.management.base import BaseCommand, CommandError
from django.db import connections

# (mbpdb_replica table, source peptide_* table, explicit shared columns --
# spelled out rather than SELECT * so this can't silently break if either
# schema's column order ever drifts)
TABLE_COPIES = [
    (
        "mbpdb_replica_proteininfo", "peptide_proteininfo",
        ["id", "header", "pid", "seq", "desc", "species"],
    ),
    (
        "mbpdb_replica_proteinvariant", "peptide_proteinvariant",
        ["id", "seq", "pvid", "protein_id"],
    ),
    (
        "mbpdb_replica_peptideinfo", "peptide_peptideinfo",
        ["id", "peptide", "protein_id", "protein_variants", "intervals", "length", "time_approved"],
    ),
    (
        "mbpdb_replica_function", "peptide_function",
        ["id", "pep_id", "function"],
    ),
    (
        "mbpdb_replica_reference", "peptide_reference",
        ["id", "func_id", "title", "authors", "abstract", "doi",
         "additional_details", "ptm", "ic50", "inhibition_type", "inhibited_microorganisms"],
    ),
]

# Delete in reverse dependency order (children before parents) to satisfy
# foreign key constraints; insert in forward order (parents before children).
DELETE_ORDER = [t[0] for t in reversed(TABLE_COPIES)]


class Command(BaseCommand):
    help = "Copy MBPDB's ProteinInfo/ProteinVariant/PeptideInfo/Function/Reference tables into the mbpdb_replica DB alias."

    def add_arguments(self, parser):
        parser.add_argument(
            "source_sqlite_path",
            help="Path to a copy of MBPDB's db.sqlite3",
        )

    def handle(self, *args, **options):
        source_path = options["source_sqlite_path"]

        connection = connections["mbpdb_replica"]
        with connection.cursor() as cursor:
            cursor.execute("ATTACH DATABASE %s AS src", [source_path])
            try:
                self.stdout.write("Clearing existing replica tables...")
                for table in DELETE_ORDER:
                    cursor.execute(f"DELETE FROM {table}")

                total = 0
                for dest_table, src_table, columns in TABLE_COPIES:
                    col_list = ", ".join(columns)
                    cursor.execute(
                        f"INSERT INTO {dest_table} ({col_list}) "
                        f"SELECT {col_list} FROM src.{src_table}"
                    )
                    total += cursor.rowcount
                    self.stdout.write(f"  {src_table} -> {dest_table}: {cursor.rowcount} rows")
            except Exception as exc:
                raise CommandError(f"Replica copy failed: {exc}")
            finally:
                cursor.execute("DETACH DATABASE src")

        self.stdout.write(self.style.SUCCESS(f"Replica loaded: {total} rows across 5 tables."))
