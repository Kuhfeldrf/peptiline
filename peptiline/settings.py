"""
Django settings for the standalone PeptiLine project.

Adapted from the MBPDB monolith's peptide/settings.py, trimmed to what
PeptiLine's three modules (data_transformation, data_analysis, heatmap_viz)
actually use (see docs/SPLIT_PLAN.md). No MBPDB models, no PEPEX/BLAST-DB
build scripts -- MBPDB search is not available in this standalone deployment
yet (see README "MBPDB integration").
"""
from dotenv import load_dotenv
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

load_dotenv(os.path.join(BASE_DIR, ".env"))

SECRET_KEY = os.environ.get("DJANGO_SECRET_KEY")

_BUILDING = os.environ.get("BUILDING", "false").lower() == "true"

if not SECRET_KEY:
    if _BUILDING:
        SECRET_KEY = "build-time-temporary-key-not-for-production"
    else:
        raise ValueError("DJANGO_SECRET_KEY environment variable is not set.")

DEBUG = os.environ.get("DEBUG", "false").lower() == "true"

_default_hosts = ["127.0.0.1", "localhost"]
_env_hosts = os.environ.get("DJANGO_ALLOWED_HOSTS", "")
_extra_hosts = [h.strip() for h in _env_hosts.split(",") if h.strip()] if _env_hosts else []
_azure_host = os.environ.get("WEBSITE_HOSTNAME", "")
if _azure_host:
    _extra_hosts.append(_azure_host)
ALLOWED_HOSTS = _default_hosts + _extra_hosts

# Celery
CELERY_WORKER_USER = "celery_user"
CELERY_BROKER_URL = "redis://127.0.0.1:6379/0"
CELERY_RESULT_BACKEND = "redis://127.0.0.1:6379/0"
CELERY_ACCEPT_CONTENT = ["json"]
CELERY_TASK_SERIALIZER = "json"
CELERY_RESULT_SERIALIZER = "json"
CELERY_TIMEZONE = "UTC"
CELERY_TASK_SOFT_TIME_LIMIT = 3600
CELERY_TASK_TIME_LIMIT = 3900

CACHES = {
    "default": {
        "BACKEND": "django_redis.cache.RedisCache",
        "LOCATION": "redis://127.0.0.1:6379/1",
        "TIMEOUT": 600,
        "OPTIONS": {"CLIENT_CLASS": "django_redis.client.DefaultClient"},
    }
}

INSTALLED_APPS = (
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    "peptiline",  # provides the `bootstrap` management command (start.sh)
    "data_transformation",
    "data_analysis",
    "heatmap_viz",
    "django_celery_progress",
    "mbpdb_replica",
)

DATABASE_ROUTERS = ["peptiline.db_router.MBPDBReplicaRouter"]

MIDDLEWARE = [
    "peptiline.middleware.HealthCheckMiddleware",
    "django.middleware.security.SecurityMiddleware",
    "whitenoise.middleware.WhiteNoiseMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

SESSION_COOKIE_AGE = 4 * 60 * 60
SESSION_SAVE_EVERY_REQUEST = True
SESSION_EXPIRE_AT_BROWSER_CLOSE = True

_default_cors = ["http://localhost:8000", "http://127.0.0.1:8000"]
_env_cors = os.environ.get("DJANGO_CORS_ORIGINS", "")
_extra_cors = [h.strip() for h in _env_cors.split(",") if h.strip()] if _env_cors else []
if _azure_host:
    _extra_cors.append(f"https://{_azure_host}")
CORS_ALLOWED_ORIGINS = _default_cors + _extra_cors

_default_csrf = ["http://localhost:8000", "http://127.0.0.1:8000"]
_env_csrf = os.environ.get("DJANGO_CSRF_ORIGINS", "")
_extra_csrf = [h.strip() for h in _env_csrf.split(",") if h.strip()] if _env_csrf else []
if _azure_host:
    _extra_csrf.append(f"https://{_azure_host}")
CSRF_TRUSTED_ORIGINS = _default_csrf + _extra_csrf

ROOT_URLCONF = "peptiline.urls"
WSGI_APPLICATION = "peptiline.wsgi.application"

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": os.path.join(BASE_DIR, "db.sqlite3"),
    },
    # Replicated MBPDB reference tables (ProteinInfo/ProteinVariant/
    # PeptideInfo/Function/Reference) -- physically separate SQLite file so
    # a `loadreplica` refresh can never lock out or corrupt PeptiLine's own
    # writes. See mbpdb_replica/ and docs/SPLIT_PLAN.md section 3.
    "mbpdb_replica": {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": os.path.join(BASE_DIR, "mbpdb_replica.sqlite3"),
    },
}

FILE_UPLOAD_HANDLERS = (
    "django.core.files.uploadhandler.TemporaryFileUploadHandler",
)
DATA_UPLOAD_MAX_MEMORY_SIZE = 104857600
FILE_UPLOAD_MAX_MEMORY_SIZE = 104857600

LANGUAGE_CODE = "en-us"
TIME_ZONE = "America/Los_Angeles"
USE_I18N = True
USE_L10N = True
USE_TZ = True
DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"

STATIC_URL = "/static/"
STATIC_ROOT = os.path.join(BASE_DIR, "static_files")
STATICFILES_DIRS = [os.path.join(BASE_DIR, "static")]
STATICFILES_STORAGE = "peptiline.storage.LenientManifestStaticFilesStorage"

MEDIA_ROOT = os.path.join(BASE_DIR, "uploads")
WORK_DIRECTORY = os.path.join(BASE_DIR, "uploads/temp")

# Headers-only protein reference used by data_transformation/heatmap_viz for
# the protein name/UniProt dictionary. MBPDB-backed BLAST search (full
# sequence FASTA + blastp) is not part of this standalone deployment yet.
PROTEIN_HEADERS_FILE = os.path.join(BASE_DIR, "protein_headers.txt")

# Species search-term translation, used by data_transformation/services/data_loader.py.
SPEC_TRANSLATE_LIST = [
    ["Human", "homo sapiens", "human", "HUMAN"],
    ["Bovine", "bos taurus", "bovine", "cow", "bovin", "BOVIN"],
    ["Sheep", "ovis aries", "sheep", "ovine", "SHEEP"],
    ["Goat", "capra hircus", "goat", "caprine", "caphi", "CAPHI"],
    ["Pig", "sus scrofa", "pig", "porcine", "sscro", "SSCRO"],
    ["Yak", "bos mutus", "yak", "mutus", "YACBA"],
    ["Rabbit", "oryctolagus cuniculus", "rabbit", "cunicu", "RABIT"],
    ["Donkey", "equus asinus", "donkey", "asinus", "EQUAS"],
    ["Camel", "camelus dromedarius", "camel", "camdr", "CAMDR"],
    ["Buffalo", "bubalus bubalis", "buffalo", "bubbu", "BUBBU"],
]

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [os.path.join(BASE_DIR, "templates")],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
            ],
        },
    },
]

SECURE_BROWSER_XSS_FILTER = True
SECURE_CONTENT_TYPE_NOSNIFF = True
X_FRAME_OPTIONS = "SAMEORIGIN"
SECURE_PROXY_SSL_HEADER = ("HTTP_X_FORWARDED_PROTO", "https")

LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "verbose": {"format": "{levelname} {asctime} {module} {message}", "style": "{"},
    },
    "handlers": {
        "console": {"class": "logging.StreamHandler", "formatter": "verbose"},
    },
    "root": {"handlers": ["console"], "level": "INFO"},
}
