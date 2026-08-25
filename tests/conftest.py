"""
pytest conftest.py — configure Django for standalone service tests.

Configures minimal Django settings via settings.configure() so `from
django.conf import settings` works without a running server, then calls
django.setup() to populate the app registry -- required since
data_transformation/services/blast_search.py imports mbpdb_replica.models
at module load, and defining a Django model class requires the app
registry to be ready (see docs/SPLIT_PLAN.md section 3 for how this
requirement got introduced).
"""
import os

os.environ.setdefault('DJANGO_SECRET_KEY', 'test-secret-key-not-for-production')

import django
from django.conf import settings

SPEC_TRANSLATE_LIST = [
    ["Bovine", "bos taurus", "btaurus", "bovin"],
    ["Human", "homo sapiens", "hsapiens", "human"],
    ["Sheep", "ovis aries", "oaries", "ovine"],
    ["Goat", "capra hircus", "chircus", "caprine"],
]

if not settings.configured:
    settings.configure(
        SECRET_KEY='test-secret-key-not-for-production',
        INSTALLED_APPS=['mbpdb_replica'],
        DATABASES={
            'default': {'ENGINE': 'django.db.backends.sqlite3', 'NAME': ':memory:'},
            'mbpdb_replica': {'ENGINE': 'django.db.backends.sqlite3', 'NAME': ':memory:'},
        },
        CACHES={'default': {'BACKEND': 'django.core.cache.backends.locmem.LocMemCache'}},
        SPEC_TRANSLATE_LIST=SPEC_TRANSLATE_LIST,
        WORK_DIRECTORY='/tmp',
        USE_TZ=True,
    )
    django.setup()
