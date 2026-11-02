# Reproduction and evidence locations

Code reuses the existing Python3.12 environment and native Julia1.9.2/JLD2 0.6.7 receipt. No new dependencies or package environments were installed. The native TSV/source hashes are revalidated, rather than decoding JLD2 with heuristics.

Run order used:

```bash
cd /home/Jack/Real_Data/lymph_node_sheaf
PY=/home/Jack/micromamba/envs/humann_py312/bin/python
"$PY" review/worm_behavior_sheaf_readiness_v1/full_cohort_qc_v5/extract_processed.py
"$PY" review/worm_behavior_sheaf_readiness_v1/full_cohort_qc_v5/audit_cohort.py
"$PY" review/worm_behavior_sheaf_readiness_v1/full_cohort_qc_v5/event_exclusions.py
"$PY" -m unittest discover -s review/worm_behavior_sheaf_readiness_v1/full_cohort_qc_v5 -p test_full_cohort.py
```

Do not rerun extraction into the populated directory: it intentionally refuses overwrites and preserves its private checkpoint. Full audit also refuses an existing primary result. For an independent reproduction use a separate project/output version or inspect/test the existing source and checksum-linked results. No source matrices are modified by any audit. Only one fish's arrays are loaded at once; the full atlas is never assembled as a dense matrix.

Private outputs are in `data/worm_behavior_sheaf_readiness_v1/full_cohort_qc_v5/`:

- `extraction_receipt_PRIVATE.json`: source archive hash, all ZIP members/CRC, extracted hashes/sizes, hashes proving original failed journal/receipts unchanged.
- `cohort_results_PRIVATE.json`: all27 animal results, split, each ROI's identity/finite/window/quantile metrics, schemas and event arrays.
- `classes_by_animal_PRIVATE.csv`, `individuals_by_animal_PRIVATE.csv`: full measured numerical/identity eligibility coverage.
- `classes_common_core_witness_PRIVATE.json`: exact witness to the joint15×15 class gate; not a sample/feature selector.
- `event_failure_decomposition_PRIVATE.json`: event-by-animal rejection counts; no thresholds altered.

Aggregate-only public tables are the extraction/clock/neural/behaviour JSON summaries, circuit/pair coverage CSVs and acceptance ledger. `TRAINING_PAIR_PANEL_COVERAGE.csv` is derived solely from training identities and numerical availability in `cohort_results_PRIVATE.json`: for each of the66 unordered predefined-circuit pairs, count intersecting training animals;≥3 permits the proposed nested animal folds. Ten pairs meet that computational eligibility rule; no performance-based selection occurred. All15/6/6 animal assignments are private and unchanged. Root current reports have an immutable previous-report snapshot in `prior_three_animal_checkpoint/`.

The primary event decomposition independently reproduces the original imported `band_events` output and verifies retained counts agree. The vectorised clock calculation uses the exact tolerances/bounds from the preserved three-animal routine. Native Julia regression tests were rerun with the existing depot/project and pass20 assertions; all Python suites pass. The longstanding frozen download fixture emits a ResourceWarning about an unclosed file, without a failed assertion. Frozen files were not changed to suppress that warning.

Readiness is not predictive validation: no transformations or models were learned, no sheaf/H0/energy constructed, no raw or individual-level files published.
