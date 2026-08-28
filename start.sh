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

# Concurrency note: this ran gunicorn's default of a single sync worker, i.e.
# exactly one request in flight at a time with everything else queued in nginx,
# while the Container App scale rule triggers at concurrentRequests=10 -- KEDA
# waited for 10 requests to stack up on a replica that could only serve one.
#
# Threads rather than extra worker processes: both databases are
# django.db.backends.sqlite3, so every additional *process* is another
# connection contending for the same write lock, and each would carry its own
# Django heap. Threads share one. CPU-heavy work goes to celery, so the GIL is
# not the binding constraint here.
#
# 8 threads also keeps /health/ answerable while a slow view is in flight,
# which the readiness probe depends on.
echo "Starting Gunicorn..."
gunicorn -b 127.0.0.1:8001 --timeout=600 \
    --worker-class gthread --workers 1 --threads 8 \
    peptiline.wsgi:application &
GUNICORN_PID=$!

# --concurrency: celery defaults to the CPU count it can see, and it reads the
# *host's* count rather than the container's cgroup limit -- it was forking 4
# prefork children, each with its own Django heap, on what was a 0.5-CPU / 1Gi
# container already peaking at ~77% of its memory limit. 2 children suits the
# 1-CPU allocation without oversubscribing the core.
#
# --max-tasks-per-child: recycle children periodically so pandas/numpy
# allocations can't accumulate across a long-lived replica. Safe for celery
# since a child is only recycled between tasks, never mid-request.
echo "Starting Celery worker..."
gosu celery_user celery -A peptiline worker --loglevel=info \
    --concurrency 2 --max-tasks-per-child 20 &
CELERY_PID=$!

wait
