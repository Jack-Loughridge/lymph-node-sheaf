# Recovery reproduction

Reuse the isolated installation; no package installation is required:

```bash
cd /home/Jack/Real_Data/lymph_node_sheaf
export JULIA_DEPOT_PATH="$PWD/data/worm_behavior_sheaf_readiness_v1/julia_tools/depot"
export JULIA_PKG_PRECOMPILE_AUTO=0
JULIA="$PWD/data/worm_behavior_sheaf_readiness_v1/julia_tools/julia-1.9.2/bin/julia"
PROJECT="$PWD/data/worm_behavior_sheaf_readiness_v1/julia_tools/project"
"$JULIA" --startup-file=no --project="$PROJECT" review/worm_behavior_sheaf_readiness_v1/export_recovery_v3/test_native_export.jl
```

The successful original native output is already preserved. To reproduce, pass a **new nonexistent** versioned output directory; the exporter rejects overwrites:

```bash
"$JULIA" --startup-file=no --project="$PROJECT" review/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2/export_identities_native.jl data/worm_behavior_sheaf_readiness_v1/recordings/NeuroPAL_labels.jld2 data/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2/native_identity_reproduction_v2
```

The reconciler currently targets the preserved first successful native output. Inspect before rerunning: it rewrites derived private summaries and its public aggregate, not source recordings. The fixed sample audit reflects three recordings; after cohort transfer use the new cohort audit first, rather than treating the sample aggregate as cohort readiness.

```bash
PY=/home/Jack/micromamba/envs/humann_py312/bin/python
"$PY" -m unittest discover -s review/worm_behavior_sheaf_readiness_v1/export_recovery_v3 -p test_recovery.py
"$PY" -m unittest discover -s review/worm_behavior_sheaf_readiness_v1 -p test_readiness.py
"$PY" -m unittest discover -s review/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2 -p test_sample_layout.py
"$PY" review/worm_behavior_sheaf_readiness_v1/export_recovery_v3/audit_transferred_cohort.py /absolute/path/to/remaining_transfer_TIMESTAMP_PRIVATE.json
```

Results this recovery:20 Julia assertions,6 new Python tests,15 frozen readiness tests and7 prior layout tests PASS. The old readiness suite emits an existing unclosed-file ResourceWarning in its download fixture; it does not alter results. Five sender/receiver SHA256+size checks and six frozen hashes PASS. No predictive model, sheaf or compatibility metric was fitted.
