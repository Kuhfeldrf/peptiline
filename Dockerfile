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

WORKDIR /app

COPY requirements.txt .
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

COPY . /app

RUN chown -R celery_user:celery_user /app && \
    chmod -R 755 /app && \
    touch /app/db.sqlite3 /app/mbpdb_replica.sqlite3 && \
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
