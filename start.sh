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

# mbpdb_replica.sqlite3 is committed to the repo already populated with MBPDB's
# reference tables, so the functional-annotation search works out of the box.
# To refresh it from a newer MBPDB dump: `manage.py loadreplica <db.sqlite3>`.

echo "Starting Gunicorn..."
gunicorn -b 127.0.0.1:8001 --timeout=600 peptiline.wsgi:application &
GUNICORN_PID=$!

echo "Starting Celery worker..."
gosu celery_user celery -A peptiline worker --loglevel=info &
CELERY_PID=$!

wait
