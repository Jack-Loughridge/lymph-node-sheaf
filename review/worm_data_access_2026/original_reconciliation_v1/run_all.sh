#!/usr/bin/env bash
set -euo pipefail
if [[ ! -f review/worm_data_access_2026/original_reconciliation_v1/reconcile.py ]]; then
  echo 'Run from /home/Jack/Real_Data/lymph_node_sheaf' >&2
  exit 1
fi
worm_reconcile_python="${WORM_QC_PYTHON:-/home/Jack/micromamba/envs/humann_py312/bin/python}"
export PYTHONDONTWRITEBYTECODE=1
"$worm_reconcile_python" review/worm_data_access_2026/original_reconciliation_v1/fetch_sources.py
"$worm_reconcile_python" -m unittest discover -s review/worm_data_access_2026/original_reconciliation_v1 -p test_reconciliation.py -v > review/worm_data_access_2026/original_reconciliation_v1/TEST_RESULTS.txt 2>&1
"$worm_reconcile_python" review/worm_data_access_2026/original_reconciliation_v1/reconcile.py "$@"
