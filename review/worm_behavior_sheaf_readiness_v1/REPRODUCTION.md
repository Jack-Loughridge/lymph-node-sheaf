# Reproduction and resumption

Run from `/home/Jack/Real_Data/lymph_node_sheaf`. The cleanup is already checkpointed: do **not** rerun the destructive entry point or overwrite its manifest. The separate host archive script is only for the exact pending read-only file, as described in MANUAL_TRANSFER.md.

No Julia executable was found; no package installation or imaging reconstruction was attempted. The existing Python environment already provides h5py/numpy. It was used read-only and not modified.

```bash
/home/Jack/micromamba/envs/humann_py312/bin/python -m unittest discover -s review/worm_behavior_sheaf_readiness_v1 -p test_readiness.py -v
/home/Jack/micromamba/envs/humann_py312/bin/python review/worm_behavior_sheaf_readiness_v1/offline_recording_audit.py
```

The second command currently records BLOCKED_NO_RECORDINGS, not a zero-neuron scientific failure. Once files exist, it records real raw/normalisation dataset shapes, numerical missingness, windows, descriptive reversal counts and source hashes privately. It does not confuse length equality with verified alignment or treat `trace_array_original` as guaranteed unprocessed measurement. Per-animal and neuron/class-by-animal tables remain private. Their current empty headers are **not measured coverage results**.

After transfer, inspect the actual `load_data_example.jl` before executing/adapting it. Use the authors' Julia `HDF5.jl`/`JLD2.jl` loading procedure to inspect actual dictionary keys/types. Do not flatten JLD2 using h5py references or interpret segmentation ROI IDs as trace-row indices. Julia/package setup, if needed, should be isolated from existing microbiome tools and logged separately. No working author-loader execution is claimed yet.

Export a private TSV from the verified Julia structure with columns:

`animal, trace_row_1based, label, neuron_class, confidence, LR`

Use tabs as delimiters. The exporter must preserve the actual individual label, class and laterality and check every trace-row index. Do not derive class names by stripping terminal letters. Document unresolved assignments rather than invent them. The dictionary schema on the public README is provisional evidence until the binary and author code are inspected.

Then:

```bash
/home/Jack/micromamba/envs/humann_py312/bin/python review/worm_behavior_sheaf_readiness_v1/offline_recording_audit.py \
 --labels-tsv data/worm_behavior_sheaf_readiness_v1/labels_PRIVATE.tsv \
 --provenance-json data/worm_behavior_sheaf_readiness_v1/verified_provenance_PRIVATE.json
```

The provenance file is a private per-recording review record, not a set of flags to assume true. Establish source signal/reconstruction masks, timestamp units, behavioural clock mapping, identity joins, circle/position units and idx_splits/session semantics from actual data and processing code. If a verified boolean observation mask exists, pass its actual HDF5 dataset path. Otherwise report the genuine-observation gate unresolved; finite values alone cannot distinguish reconstruction.

Reversal table counts are descriptive only. After clock and boundary calibration are verified, derive copper-entry event counts and temporally censored behavioural targets using the frozen definition in VALIDATION_DESIGN.md; do not guess a conversion from the README's example. Confirm variation and repeated-event animals before running a diagnostic.

The 27-animal split was frozen from public filenames and is stored under ignored data/. Publish only group counts/hash. Do not acquire the remaining 24 processed HDF5 files until the three-file audit and dictionary join pass. The current directory intentionally contains no model-fitting script: biological prerequisites have not passed. Resume the next scientific step after the small transfer; no broad search is needed.
