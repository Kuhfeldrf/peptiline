"""
pytest conftest.py — configure Django for standalone service tests.

Strategy:
  - Stub out only `peptide.celery` in sys.modules BEFORE the real `peptide`
    package is imported.  When peptide/__init__.py does
    `from .celery import app as celery_app` Python finds the stub first,
    so no Redis connection is attempted.
  - Configure Django settings via settings.configure() so `from django.conf
    import settings` works without a running server.
"""
import os
import sys
import types
from unittest.mock import MagicMock

# ---------------------------------------------------------------------------
# 1. Stub out peptide.celery BEFORE any package import
#    This prevents the real celery.py from trying to connect to Redis.
# ---------------------------------------------------------------------------
_celery_stub = types.ModuleType('peptide.celery')
_celery_stub.app = MagicMock()
sys.modules['peptide.celery'] = _celery_stub

# ---------------------------------------------------------------------------
# 2. Configure Django with minimal in-memory settings
# ---------------------------------------------------------------------------
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
        INSTALLED_APPS=[],
        DATABASES={},
        CACHES={'default': {'BACKEND': 'django.core.cache.backends.locmem.LocMemCache'}},
        SPEC_TRANSLATE_LIST=SPEC_TRANSLATE_LIST,
        WORK_DIRECTORY='/tmp',
        USE_TZ=True,
    )
