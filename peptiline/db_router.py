class MBPDBReplicaRouter:
    """
    Pins every mbpdb_replica model to the `mbpdb_replica` DB alias, and
    everything else to `default`. Keeps the replicated MBPDB tables in a
    physically separate SQLite file from PeptiLine's own operational data
    (sessions, admin, django_celery_progress) -- see docs/SPLIT_PLAN.md
    section 3, "replica tables should live in a separate SQLite file".
    """

    replica_app = "mbpdb_replica"

    def db_for_read(self, model, **hints):
        return self.replica_app if model._meta.app_label == self.replica_app else None

    def db_for_write(self, model, **hints):
        return self.replica_app if model._meta.app_label == self.replica_app else None

    def allow_relation(self, obj1, obj2, **hints):
        replica_labels = {obj1._meta.app_label == self.replica_app, obj2._meta.app_label == self.replica_app}
        if len(replica_labels) == 1:
            return True
        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        if app_label == self.replica_app:
            return db == self.replica_app
        return db != self.replica_app
