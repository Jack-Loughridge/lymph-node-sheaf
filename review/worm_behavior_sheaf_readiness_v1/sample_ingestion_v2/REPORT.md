# Candidate C: verified transfer and actual three-recording audit

**Transfer integrity PASS: all five Windows SHA-256 hashes and byte counts match the VM files.** The original files and receipt remain private and unchanged. Actual numerical sample QC passes the tests performed. Full technical/biological readiness is **PARTIAL GO / pending native identity and observation-provenance checks**, not a full-model GO.

## What was actually inspected

The five transferred payloads were opened/read, including the actual author Julia loading code. There are three HDF5 recordings, each with 30 numerical datasets and 1,600 neural time points. Across them are 439 ROI-recordings, 139–150 per animal, **not 439 independent neurons or animals**. Their durations are 961.896–961.918 seconds. Three nominal independent animals are represented; this does not establish the required 15-animal cohort.

All original, z-scored, F/Fmean and F/F20 arrays are finite, and no original trace is constant. Across 439 traces there are 698,010 valid **numerical** ten-history forecast windows. These are not independent biological replicates or certified genuine-observation windows. Native neural intervals range 0.4475–0.7250 seconds, with median approximately 0.600019 seconds; no interval exceeds 1.5 median intervals. There are 1,220 NIR intervals exceeding the analogous frozen threshold; these remain in the timing audit rather than being silently interpolated away.

The actual storage is **time × neurons in h5py**, while the author Julia example uses neurons × time. Event and segment matrices are likewise reversed. A new, tested adapter corrects that implementation convention; it does not change eligibility, confidence, masking, split or performance thresholds. Frozen v1 code/protocol files remain byte-identical.

The exported z-scores match per-neuron full-record sample-standard-deviation scaling (ddof=1), and F/Fmean and F/F20 match full-record mean and twentieth-percentile division. These are recording-wide derived transformations, not training-only normalisation. No model has used them. The `_original` name does not prove absence of earlier processing.

## Timing and behavioural evidence

Velocity, reversals, head angle, angular velocity, pumping and supplied distance-to-copper-centre vectors are all finite in all three recordings. There are 32–42 stored reversal intervals per recording, 113 total, with nonconstant velocity and other continuous behavioural variables. Two recordings begin during a reversal; the initial intervals are left-censored. There are 111 new reversal onsets after the recording starts. These are reversal counts, **not verified copper-encounter counts**.

The actual table indices reproduce the inspected author exporter's `diff` convention, including its special initial-reversal case. A generic 'subtract one from every event index' rule would misinterpret these boundaries. Canonical reversal states should be derived from the binary vector with censoring preserved.

An explicit metadata reconciliation found contiguous NIR bins whose mean timestamp and mean x/y stage positions jointly reproduce the corresponding stored confocal-frame fields. All 4,800 frames match; 4,799 have a unique matching bin and one is ambiguous. This is strong exported-clock alignment evidence, not mere equality of array lengths. Original DAQ triggers and the original confocal-to-NIR index vector are not exported, so the result is a derived reconstruction and not a claim of recovered ground-truth hardware intervals. Ambiguity must be masked or retained as timing uncertainty. Do not use binned behaviour as an instantaneous predictor before its acquisition bin has ended.

The public source defines 10,000 stage units per millimetre, but restoring the exported stage-origin offsets alone does **not** reproduce supplied radial distances from circle centre. Head/camera coordinates or exact candidate-specific geometry code remain needed to independently reconstruct that distance. No unrecorded coordinate correction was fitted or silently assumed, and no boundary-event or avoidance forecast was run.

## Genuine measurements versus processed values

The public ANTSUN notebook explicitly invokes `process_traces` with interpolation enabled. The CaAnalysis implementation can interpolate activity and marker traces, divide channels and fit bleaching corrections; the historical HDF5 exporter names processed `traces_array` as `trace_array_original`. The actual files lack per-neuron/per-time detection or reconstruction masks. The uninspected native identity metadata may expose additional provenance; irrecoverability is not established. Candidate-specific export additions differ from the historical exporter, so the precise pipeline/version used for these files is **not established** solely by that code.

Thus all-finite signals are numerically usable but do not establish the genuine-observation fraction. No exact-linear interior points were found at machine-precision tolerance; that cannot rule out interpolation before nonlinear ratio/bleach correction. No such detector was used as an observation mask.

The smallest useful provenance request is an observed/detected mask for each final trace row in both activity and reference channels, the ROI-to-final-row mapping, preprocessing/version/bleach parameters, and, if necessary, pre-interpolation channel traces. This is a request for small processed exports, not raw microscopy reconstruction. Recording-wide future-dependent correction also needs review before claiming causal short-history forecasting.

## Native identity dependency and installation

The original JLD2 opens as an HDF5 container and has the author-named `dict_match_dict` root. **Its Julia objects have not been decoded heuristically in Python.** The actual author loader uses `jldopen`, indexes that root, closes the file, then queries the two row/class lookup views. No installed Julia runtime or Python Julia bridge was found. An official-runtime request failed DNS. This executor also lacks the user's configured `tda-vm` SSH alias, so normal-SSH network reachability could not be tested from it.

[JULIA_SETUP.md](JULIA_SETUP.md) provides the normal SSH command, official Linux x86_64 archive/checksum, Windows/SCP fallback and isolated package/depot setup. Only JLD2 is installed; HDF5.jl, imaging/GPU dependencies and changes to existing Python environments are unnecessary. The native exporter follows the actual author loading operation and cross-checks both redundant identity views before writing private labels and a version/hash receipt. The installer is bash-syntax checked, but neither Julia execution nor package installation is claimed completed.

Do not substitute `neuropal_registration/roi_match_confidence` for manual label confidence: its actual numeric range is not the dictionary's 1–5 scale. The original/skip trace mappings are identity relations stored as unordered pairs; they must not be treated as a new ordering of trace rows. Individual left/right assignments and class overlap remain unverified until the native dictionary is interpreted.

## Decision and next steps

1. Run the isolated Julia installer through the normal `tda-vm` alias, or transfer the verified official Linux archive using the documented fallback. The five dataset files do not need downloading again.
2. Run native JLD2 validation, then reconcile identities against each recording's actual trace-row bounds and preserve uncertainty/conflicts. Quantify common classes and side-resolved identities without changing the >=15-animal/>=15-class gates.
3. Resolve or explicitly bound original-observation/preprocessing and boundary-geometry provenance. Numerical integrity alone is not the full sample technical gate.
4. Only if that sample gate passes, acquire the remaining processed recordings and perform the frozen readiness/diagnostic protocol. The 15/6/6 split and evaluation rules remain unchanged; all older results are preserved.

**No full-collection acquisition, neural model, sheaf, compatibility metric or behavioural forecast was run.** The study has not received a scientific NO-GO. The remaining blocks are specific, actionable dependencies, not generic publication-metadata uncertainty.

To reproduce the numerical checks from the project root:

```bash
/home/Jack/micromamba/envs/humann_py312/bin/python review/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2/reconcile_clock_bins.py
/home/Jack/micromamba/envs/humann_py312/bin/python review/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2/sample_numerical_audit.py
```

Seven new layout/index tests pass; the existing fifteen tests still pass. Native Julia validation remains pending and is not counted as a passed test. All actual animal/neuron tables, array schemas, receipt details and source-file hashes remain under ignored data/. Public outputs are code and aggregate research findings only.
