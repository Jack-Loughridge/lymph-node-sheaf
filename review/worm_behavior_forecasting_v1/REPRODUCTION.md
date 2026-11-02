# Reproduce without opening validation or final-test recordings

Use the existing Python environment; no installation was required:

```bash
cd /home/Jack/Real_Data/lymph_node_sheaf
PY=/home/Jack/micromamba/envs/humann_py312/bin/python
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
"$PY" -m unittest discover -s review/worm_behavior_forecasting_v1 -p 'test_*.py'
"$PY" review/worm_behavior_forecasting_v1/run_diagnostic.py
"$PY" review/worm_behavior_forecasting_v1/summarise.py
"$PY" review/worm_behavior_forecasting_v1/audit_predictions.py
"$PY" review/worm_behavior_forecasting_v1/make_report.py
MPLCONFIGDIR=/tmp/worm_forecast_matplotlib "$PY" review/worm_behavior_forecasting_v1/plot_aggregates.py
```

The runner resumes only checksum-matching pair/history/group checkpoints; changed core code/protocol requires a versioned run, never overwriting old fitted results. Complete checkpoints do not refit. Pre-fit implementation hashes are immutable and independently audited. No observations, preprocessing parameters or performance from the six validation/six final-test recordings are loaded. The only all-cohort inputs are the pre-existing split mapping and native identity metadata, filtered to training before recording access; the public training-only support table is checked exactly.

Private inputs: `data/worm_behavior_sheaf_readiness_v1/recordings/`, the existing native Julia TSV/receipt, private split mapping; independent template `data/bioelectric_focused_feasibility_2026/white_whole.tsv`. Native JLD2 is not heuristically decoded or reinstalled. The corrected Julia export was already validated; its TSV SHA256 is enforced here. Source recordings are read-only. The source list in `run_receipt_PRIVATE.json` pins SHA256 for every training HDF5 opened.

Private outputs: `data/worm_behavior_forecasting_v1/` contains the run receipt,60 panel/history/group checkpoints (including three unevaluable day entries),261 prediction NPZs, complete nested fold metrics and the animal-level table. Checkpoints record fitting animal IDs, chosen alphas, inner scores, apparent training errors, outer errors and prediction-file hashes. Raw source data are not copied. Fitted coefficients/scalers are deterministically recoverable from the recorded training sets, selected alphas, source hashes and frozen code; they are not published. All10 eligible pair panels and both directions are reported. No file under this private output directory is on the GitHub publication allowlist.

Public outputs:

- `MODEL_COMPARISONS.csv`:798 aggregate model/direction/history/group rows, original-unit/normalised and apparent-training/outer MSE.
- `DIRECTED_CONTRASTS.csv`:1,824 complete aggregate contrasts across both scales and all frozen sensitivities; animal and date-cluster intervals, support and influence.
- `ALL_DIRECTIONS_PRIMARY.csv`:all20 directions at ten-sample history, no favourable selection.
- `POOLED_CONTRASTS.csv`:equal-animal pooled estimates; `DAY_MATCHED_POOL.csv` separates date exclusion from loss of unsupported panels.
- `ELIGIBILITY_ALL_66.csv`, `PAIR_EVALUABILITY.csv`:unsupported pair/date states remain explicit.
- `ANATOMICAL_COMPARATOR.csv`:reference edge masks by generic neuron labels, with no individual-animal IDs.
- `IMPLEMENTATION_PREREGISTRATION.md`, `PREREGISTRATION_RECEIPT.json`, `PREFIT_CODE_RECEIPT.json`:pre-fit constraints and code fingerprints.
- `FORECAST_DIAGNOSTIC.png/.pdf`:aggregate-only exportable figure; no animal-level traces or IDs.

No cell/time-window pseudo-replication is used. The ordinary and day sensitivity predictions can overlap; they are never counted as new independent animals. Bootstrap resamples the observed cross-validation losses, preserving all panels for an animal; it does not refit the algorithm inside each bootstrap and therefore is a descriptive conditional uncertainty analysis. Per-edge intervals are unadjusted; no causal/discovery claims follow.

The original six-file frozen criteria/specification/code hashes and all source files remain unchanged. Publication uses an explicit allowlist restricted to this review directory, excluding all raw/animal data. The historical readiness reports and failed-transfer records remain intact.
