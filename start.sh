#!/bin/bash
set -e

cleanup() {
    echo "Cleaning up processes..."
    kill $GUNICORN_PID 2>/dev/null || true
    kill $CELERY_PID 2>/dev/null || true
    service nginx stop || echo "Failed to stop Nginx"
    exit 0
}

trap cleanup SIGTERM SIGINT

echo "Starting Redis server..."
service redis-server start || echo "Redis server failed to start"

echo "Starting Nginx..."
nginx -t && service nginx start || { echo "Nginx failed to start: $(nginx -t 2>&1)"; exit 1; }

cd /app

echo "Running startup bootstrap (migrate, clearsessions, superuser check)..."
python manage.py bootstrap || echo "WARNING: bootstrap step failed"

# Load the bundled MBPDB data snapshot into mbpdb_replica so search works
# out of the box. mbpdb_seed.sqlite3 is a checked-in copy, refreshed by
# rebuilding the image from a newer MBPDB db.sqlite3 -- see
# docs/SPLIT_PLAN.md section 3 for the real scheduled-refresh transport
# this is standing in for.
echo "Loading MBPDB replica data..."
python manage.py loadreplica /app/mbpdb_seed.sqlite3 || echo "WARNING: loadreplica failed"

echo "Starting Gunicorn..."
gunicorn -b 127.0.0.1:8001 --timeout=600 peptiline.wsgi:application &
GUNICORN_PID=$!

echo "Starting Celery worker..."
gosu celery_user celery -A peptiline worker --loglevel=info &
CELERY_PID=$!

wait
