#!/usr/bin/env bash
# Build OptimOTU_v6.sif (default). Uses GITHUB_PAT from the
# environment so the token is only written into a temp file during build and
# is never stored in the image.
# Run from the project root:
#   ./apptainer/build.sh          # v6
#   ./apptainer/build.sh v[n]
# Optional: set OPTIMOTU_RENV_CACHE to override the default host renv cache
# path ($HOME/.cache/R/renv). When the directory exists, it is bind-mounted to
# /renv/cache for build-time reuse and is not copied into the image.
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

VERSION="${1:-v6}"
DEF_FILE="apptainer/OptimOTU_${VERSION}.def"
SIF_FILE="apptainer/OptimOTU_${VERSION}.sif"
if [ ! -f "$DEF_FILE" ]; then
  echo "Container definition '$DEF_FILE' not found."
  exit 1
fi

TOKEN_FILE="$PROJECT_ROOT/apptainer/GITHUB_TOKEN"
if [ -n "${GITHUB_PAT:-}" ]; then
    printf '%s' "$GITHUB_PAT" > "$TOKEN_FILE"
    echo "Using GITHUB_PAT for build (token file created)."
else
    touch "$TOKEN_FILE"
    echo "GITHUB_PAT not set; creating empty token file (you may hit GitHub rate limits)."
fi

cleanup() { rm -f "$TOKEN_FILE"; }
trap cleanup EXIT

echo "Building ${SIF_FILE} from ${DEF_FILE}"
RENV_CACHE_DIR="${OPTIMOTU_RENV_CACHE:-$HOME/.cache/R/renv}"
BIND_ARGS=()
if [ -d "$RENV_CACHE_DIR" ]; then
    BIND_ARGS+=(--bind "${RENV_CACHE_DIR}:/renv/cache")
    echo "Using host renv cache bind: ${RENV_CACHE_DIR} -> /renv/cache"
else
    echo "Host renv cache not found at ${RENV_CACHE_DIR}; building without bind."
fi

apptainer build "${BIND_ARGS[@]}" \
    "$SIF_FILE" \
    "$DEF_FILE"
