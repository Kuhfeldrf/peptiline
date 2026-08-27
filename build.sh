#!/bin/bash

set -euo pipefail

IMAGE_NAME=peptiline
CONTAINER_NAME=peptiline_dev

# Usage help
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  echo "Usage: $0 [clean] [--mount]"
  echo "  clean   Perform a full rebuild without cache and prune old container/image"
  echo "  --mount Bind-mount source into container for rapid iteration (dev only)"
  exit 0
fi

CLEAN_BUILD=false
BIND_MOUNT=false
for arg in "$@"; do
  case "$arg" in
    clean) CLEAN_BUILD=true ;;
    --mount) BIND_MOUNT=true ;;
  esac
done

# Stop/remove only our named container if it exists
sudo docker rm -f "$CONTAINER_NAME" 2>/dev/null || true

# Build with cache by default; use --no-cache for clean builds
if $CLEAN_BUILD; then
  echo "Performing clean build (no cache)..."
  DOCKER_BUILDKIT=1 sudo docker build --no-cache -t "$IMAGE_NAME" .
else
  echo "Performing incremental build (using cache)..."
  DOCKER_BUILDKIT=1 sudo docker build -t "$IMAGE_NAME" .
fi

# Run container; expose only the Django/Nginx front port required externally
RUN_ARGS=(-p 8000:8000 --name "$CONTAINER_NAME")

if $BIND_MOUNT; then
  echo "Starting with bind mounts for rapid iteration..."
  RUN_ARGS+=(
    -v "$(pwd):/app:rw"
  )
fi

sudo docker run "${RUN_ARGS[@]}" "$IMAGE_NAME"
