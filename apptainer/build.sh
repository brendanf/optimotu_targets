#!/usr/bin/env bash
# Build OptimOTU_v7.sif (default) or OptimOTU_v6.sif. Uses GITHUB_PAT from the
# environment so the token is only written into a temp file during build and
# is never stored in the image.
#
# GITHUB_PAT is required: renv resolves GitHub packages via api.github.com,
# and anonymous requests are limited to ~60/hour (easy to exhaust on rebuilds).
# Create a classic PAT with public_repo (or fine-grained Contents: Read on the
# needed repos) and export it before building.
#
# Run from the project root:
#   export GITHUB_PAT=ghp_...
#   ./apptainer/build.sh          # v7
#   ./apptainer/build.sh v7
#   ./apptainer/build.sh v6
# Optional: set OPTIMOTU_RENV_CACHE to override the default host renv cache
# path ($HOME/.cache/R/renv). When the directory exists, it is bind-mounted to
# /renv/cache for build-time reuse and is not copied into the image.
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

VERSION="${1:-v7}"
DEF_FILE="apptainer/OptimOTU_${VERSION}.def"
SIF_FILE="apptainer/OptimOTU_${VERSION}.sif"
if [ ! -f "$DEF_FILE" ]; then
  echo "Container definition '$DEF_FILE' not found."
  exit 1
fi

TOKEN_FILE="$PROJECT_ROOT/apptainer/GITHUB_TOKEN"
if [ -z "${GITHUB_PAT:-}" ]; then
    echo "ERROR: GITHUB_PAT is not set." >&2
    echo "renv needs authenticated GitHub API access to install packages" >&2
    echo "(e.g. alessandrozito/BayesANT). Export a PAT and retry:" >&2
    echo "  export GITHUB_PAT=ghp_..." >&2
    exit 1
fi
printf '%s' "$GITHUB_PAT" > "$TOKEN_FILE"
echo "Using GITHUB_PAT for build (token file created)."

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
