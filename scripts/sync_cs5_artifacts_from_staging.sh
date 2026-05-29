#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DEFAULT_STAGING_DATA="/Users/othomas/Desktop/CASE STUDY 5 SWEEP/case_study_artifacts/data"
CS5_STAGING_DATA="${CS5_STAGING_DATA:-${DEFAULT_STAGING_DATA}}"

if [[ ! -d "${CS5_STAGING_DATA}" ]]; then
  echo "Could not find staged artifact data directory:" >&2
  echo "  ${CS5_STAGING_DATA}" >&2
  echo >&2
  echo "Set CS5_STAGING_DATA to the directory containing artifact_manifest.json." >&2
  exit 1
fi

if [[ ! -f "${CS5_STAGING_DATA}/artifact_manifest.json" ]]; then
  echo "Staged artifact directory is missing artifact_manifest.json:" >&2
  echo "  ${CS5_STAGING_DATA}" >&2
  exit 1
fi

rsync -a --delete --exclude ".DS_Store" "${CS5_STAGING_DATA}/" "${REPO_ROOT}/data/"

echo "Synced Case Study 5 staged artifacts into:"
echo "  ${REPO_ROOT}/data"
