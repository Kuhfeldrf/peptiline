"""
Read-only replica of the four MBPDB tables PeptiLine's BLAST search needs
(see docs/SPLIT_PLAN.md section 3). Field-for-field copy of MBPDB's
peptide.models schema for ProteinInfo/ProteinVariant/PeptideInfo/Function/
Reference -- kept in this separate app so a Django database router can pin
all of it to the `mbpdb_replica` DB alias, physically apart from PeptiLine's
own operational SQLite file. MBPDB is the sole writer; nothing in PeptiLine
ever creates/updates rows here, only `loadreplica` (via `loaddata`) and
reads from data_transformation/services/blast_search.py.
"""
from django.db import models


class ProteinInfo(models.Model):
    header = models.CharField(max_length=1000)
    pid = models.CharField(max_length=30, db_index=True)
    seq = models.CharField(max_length=10000)
    desc = models.CharField(max_length=500)
    species = models.CharField(max_length=150)


class ProteinVariant(models.Model):
    seq = models.CharField(max_length=10000)
    pvid = models.CharField(max_length=30)
    protein = models.ForeignKey(ProteinInfo, on_delete=models.CASCADE, related_name="orig_proteins")


class PeptideInfo(models.Model):
    # Indexed: blast_search.py filters on `peptide` once per exact-match sequence
    # in an uploaded dataset (thousands of lookups per annotation run).
    peptide = models.CharField(max_length=500, db_index=True)
    protein = models.ForeignKey(ProteinInfo, on_delete=models.CASCADE, related_name="proteins")
    protein_variants = models.CharField(max_length=100, default='')
    intervals = models.CharField(max_length=100)
    length = models.IntegerField()
    time_approved = models.DateTimeField()


class Function(models.Model):
    pep = models.ForeignKey(PeptideInfo, related_name="functions", on_delete=models.CASCADE)
    function = models.CharField(max_length=400)

    class Meta:
        unique_together = [['pep', 'function']]


class Reference(models.Model):
    func = models.ForeignKey(Function, related_name="references", on_delete=models.CASCADE)
    title = models.CharField(max_length=300)
    authors = models.CharField(max_length=300)
    abstract = models.CharField(max_length=1000)
    doi = models.CharField(max_length=100)
    additional_details = models.CharField(max_length=400, default='')
    ptm = models.CharField(max_length=200, default='')
    ic50 = models.FloatField(null=True, blank=True)
    inhibition_type = models.TextField(null=True, blank=True)
    inhibited_microorganisms = models.TextField(null=True, blank=True)
