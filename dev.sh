#!/usr/bin/env bash
# dev.sh — start the Django dev server for local testing
# Usage:  bash dev.sh [port]   (default port: 8000)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${1:-8000}"

cd "$REPO_ROOT"

# ── Env vars required by settings.py ─────────────────────────────
# .env (git-ignored) is loaded automatically by settings.py via python-dotenv;
# these exports are just a fallback so the script also works without one.
export DJANGO_SECRET_KEY="${DJANGO_SECRET_KEY:-dev-secret-key-not-for-production}"
export DJANGO_ALLOWED_HOSTS="${DJANGO_ALLOWED_HOSTS:-localhost,127.0.0.1}"
export DEBUG="${DEBUG:-True}"

# ── Sanity checks ─────────────────────────────────────────────────
echo "=== Django system check ==="
python3 manage.py check --deploy 2>/dev/null || python3 manage.py check

echo ""
echo "=== Applying migrations ==="
python3 manage.py migrate --run-syncdb --noinput
python3 manage.py migrate --database=mbpdb_replica --noinput

# mbpdb_replica lives in its own sqlite file (settings.py DATABASES alias) and
# starts empty until loaded — without this, MBPDB search silently returns 0
# results instead of erroring, which looks like a bug rather than empty data.
# Reloaded unconditionally on every start (same as start.sh in the container)
# since loadreplica just replaces the tables wholesale — cheap and idempotent.
if [ -f mbpdb_seed.sqlite3 ]; then
  echo ""
  echo "=== Loading MBPDB replica seed data ==="
  python3 manage.py loadreplica mbpdb_seed.sqlite3
fi

echo ""
echo "=== Collecting static files ==="
python3 manage.py collectstatic --noinput --clear 2>&1 | tail -3

echo ""
echo "=== Starting dev server on http://127.0.0.1:$PORT ==="
echo ""
echo "  Data Analysis:       http://127.0.0.1:$PORT/data_analysis/"
echo "  Data Transformation: http://127.0.0.1:$PORT/data_transformation/"
echo "  Heatmap:             http://127.0.0.1:$PORT/heatmap/"
echo ""
echo "  Press Ctrl+C to stop."
echo ""

exec python3 manage.py runserver "127.0.0.1:$PORT"
