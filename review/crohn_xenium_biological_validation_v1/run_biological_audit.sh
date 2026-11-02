#!/usr/bin/env bash
set -euo pipefail
if [[ ! -f review/crohn_xenium_biological_validation_v1/AUDIT_PROTOCOL.json ]]; then
  echo 'Run from /home/Jack/Real_Data/lymph_node_sheaf' >&2
  exit 1
fi
crohn_audit_python="${CROHN_AUDIT_PYTHON:-/home/Jack/micromamba/envs/humann_py312/bin/python}"
export PYTHONDONTWRITEBYTECODE=1
python3 review/crohn_xenium_biological_validation_v1/extract_workbook.py
"$crohn_audit_python" -m unittest discover -s review/crohn_xenium_biological_validation_v1 -p test_spatial_audit.py -v > review/crohn_xenium_biological_validation_v1/SPATIAL_TEST_RESULTS.txt 2>&1
"$crohn_audit_python" review/crohn_xenium_biological_validation_v1/spatial_audit.py > data/crohn_xenium/biological_validation_v1/spatial_audit.log
"$crohn_audit_python" review/crohn_xenium_biological_validation_v1/coverage_sensitivity.py
"$crohn_audit_python" review/crohn_xenium_biological_validation_v1/liberal_bound_sensitivity.py
"$crohn_audit_python" review/crohn_xenium_biological_validation_v1/subtype_upper_bounds.py
"$crohn_audit_python" review/crohn_xenium_biological_validation_v1/additional_marker_audit.py > data/crohn_xenium/biological_validation_v1/additional_marker_audit.log
