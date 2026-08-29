# Dockerfile
FROM python:3.10

ENV PYTHONUNBUFFERED=1
ENV DJANGO_SETTINGS_MODULE=peptiline.settings
ENV PIP_ROOT_USER_ACTION=ignore

# System dependencies + gosu (privilege dropping) + celery_user, mirroring
# the MBPDB monolith's Dockerfile. ncbi-blast+ is needed here now that
# blast_search.py queries the mbpdb_replica database directly (see
# docs/SPLIT_PLAN.md section 3, database replication) -- no perl/PEPEX
# scripts, those stay MBPDB-only.
RUN apt-get update && apt-get install -y \
    gosu \
    nginx \
    dos2unix \
    nano \
    sqlite3 \
    ncbi-blast+ \
    redis-server \
    build-essential \
    curl \
    && useradd -r -s /sbin/nologin celery_user \
    && rm -rf /var/lib/apt/lists/*

# Shared libraries headless Chrome needs at runtime (see kaleido_get_chrome
# below) — python:3.10 doesn't ship these by default.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libnspr4 \
    libnss3 \
    libdbus-1-3 \
    libatk1.0-0 \
    libatk-bridge2.0-0 \
    libcups2 \
    libxkbcommon0 \
    libatspi2.0-0 \
    libxcomposite1 \
    libxdamage1 \
    libxfixes3 \
    libxrandr2 \
    libgbm1 \
    libasound2 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Downloads a self-contained Chrome for Testing build that Kaleido locates
# automatically at runtime (no system Chrome/PATH entry needed) — required
# for the PNG/SVG plot export endpoints in data_analysis and heatmap_viz.
RUN kaleido_get_chrome -f

# matplotlib's first import scans every system font to build its cache
# (fontManager) if one isn't already on disk -- a multi-second cost that,
# without this, was paid by whichever request first hit a heatmap render
# after every scale-from-zero cold start. Baking the cache into the image
# here means it's just read, not rebuilt, at request time. Runs as root,
# same as gunicorn does at runtime, so the cache lands where it's looked for.
RUN python -c "import matplotlib.pyplot"

COPY . /app

# mbpdb_replica.sqlite3 is the checked-in, populated snapshot the BLAST /
# functional-annotation search reads -- it MUST arrive via `COPY . /app` above
# (nothing loads it at deploy time). Fail the build loudly if it is missing or
# was truncated (e.g. re-added to .dockerignore), rather than shipping an image
# whose every database search fails at runtime. db.sqlite3, by contrast, is
# created fresh by `migrate` on first boot.
RUN test -s /app/mbpdb_replica.sqlite3 \
    && [ "$(stat -c%s /app/mbpdb_replica.sqlite3)" -gt 100000 ] \
    || (echo 'ERROR: mbpdb_replica.sqlite3 missing/empty in build context -- check .dockerignore' && exit 1)

RUN chown -R celery_user:celery_user /app && \
    chmod -R 755 /app && \
    touch /app/db.sqlite3 && \
    chown celery_user:celery_user /app/db.sqlite3 /app/mbpdb_replica.sqlite3 && \
    chmod 664 /app/db.sqlite3 /app/mbpdb_replica.sqlite3 && \
    mkdir -p /app/uploads/temp && \
    chown celery_user:celery_user /app/uploads/temp && \
    chmod 750 /app/uploads/temp

COPY start.sh /app/start.sh
RUN chmod +x /app/start.sh

COPY nginx.conf /etc/nginx/nginx.conf

# Collect static files (BUILDING=true allows this without SECRET_KEY)
RUN BUILDING=true python manage.py collectstatic --noinput --verbosity 1 && \
    chmod -R 755 /app/static_files

EXPOSE 8000

CMD ["/app/start.sh"]
