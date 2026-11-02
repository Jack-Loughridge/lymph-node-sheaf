# Candidate C: storage cleanup and dynamical-sheaf readiness

**Current decision: no modelling GO. Readiness is NOT EVALUABLE because no new physiological recording could be acquired.** At most PARTIAL GO applies to acquisition/preparation; a smaller usable circuit has not yet been demonstrated. There is no evidence-backed scientific NO-GO either. Network failure is not dataset inadequacy. The workflow remains incomplete at the representative-data gate.

## 1. Cleanup completed within the executor's writable scope

| Exact source download | Verified bytes | Outcome |
|---|---:|---|
| `/home/Jack/Real_Data/lymph_node_sheaf/data/crohn_xenium/scp3259/xenium5k_scVI_CellCharter_batch2_reannot_SCP.h5ad` | 8,984,620,838 | **Removed** after size, SHA-256, path-component and open-handle checks |
| `/home/Jack/Real_Data/spatial_ibd_sheaf/archives/GSE312415_files.tar.gz` | 6,464,275,219 | **Preserved:** unlink failed with errno 30, read-only filesystem |

The removed H5AD's SHA-256 matched the frozen value `74759d2d474964dcb537533127f0683bdb9180030f3b94a8f3f8271a061494bb`. The archive's newly calculated SHA-256 is `32fa2857563a117f974ed507083a996c34a260b52633b30fc3a9e1bc334c7cc6`. Its retained member inventory and existing extraction code establish that it is a downloaded GEO payload bundle; no archive extraction was performed. Restoring an equivalent bundle may not recreate its exact compression bytes.

**Removed logical data: 8.985 GB (8.368 GiB).** Allocated blocks removed: 8,985,174,016 bytes. Observed increase in free space at the cleanup checkpoint: 8,984,469,504 bytes. Free space then: **18,255,200,256 bytes (18.26 GB / 17.00 GiB)**. Minor report/cache writes subsequently reduce this slightly. This is sufficient for the approximately 320 MB processed worm collection; the pending archive could reclaim a further 6.46 GB.

`cleanup_manifest.json` is the durable, per-file before/after journal with exact paths, recoverability, public source URLs and checksums. Only two literal allowlisted paths were eligible. Every ancestor was checked for symlinks; only regular, singly linked files were permitted. Neither broad recursive deletion nor privilege/network bypass was used. Proc/lsof checks found no accessible users of either target. Process inspection is limited to executor-visible processes, with lsof filesystem warnings recorded; no host-wide guarantee is claimed.

Largest top-level directories before cleanup included `humann_databases` (79.16 GB), `AGORA2_SBML_clean` (50.88 GB), `external_validation` (12.90 GB) and `out_graphs` (9.79 GB). These belong to other/possibly active work and were **not deletion targets**. Their complete size inventory is private. The 323.9 MB Crohn gene-subset cache is a derived research artifact and was preserved. Other uncertain datasets were left alone, not queued for deletion without approval.

All code, repositories, specifications, tests, reference dictionaries and aggregate research results remain. Fifty-three frozen prior review/specification hashes are unchanged. The three existing worm NWBs and patient-linkage workbook are present and privately hashed. Original-export work remains intact. The microbiome PhD project and source data were not modified. Historical reports describe the earlier audit state; this cleanup manifest now documents the intentional source-download removal and how it can be restored.

The exact pending archive can be removed from a normal `tda-vm` session using the guarded script in [MANUAL_TRANSFER.md](MANUAL_TRANSFER.md). Its failure here was a filesystem boundary, not an automatic approval rejection. The manifest is not overwritten by that optional host-side action; a separate receipt records it.

## 2. Acquisition: real manifest, blocked binary payloads

The [public Dryad file listing](https://datadryad.org/dataset/doi:10.5061/dryad.w9ghx3g4v) contains 27 named per-animal HDF5 files plus a shared NeuroPAL dictionary and author Julia loader. Listed total: 320.03 MB; individual recordings 9.40–12.45 MB, dictionary 41.61 MB. The three earliest chronological files were fixed as the sample before outcomes were seen. PUBLIC_FILE_MANIFEST.csv records names and rounded sizes; repository integrity checksums were not available in the accessible listing.

After the permitted cleanup finished, one direct anonymous sample request failed DNS in this executor. The remaining same-host shell requests were not retried. The web reader returned 403 for the manifest-linked sample, dictionary and loader endpoints. The runtime environment-status tool/network-policy snapshot were unavailable here, so this establishes neither VM-wide network failure nor a permanent public access restriction. No file payload or raw video was downloaded. The complete collection was not requested because the sample gate did not pass.

[MANUAL_TRANSFER.md](MANUAL_TRANSFER.md) provides exact manifest-linked download URLs and Windows instructions. `TRANSFER_WINDOWS.ps1` uses **`tda-vm`**, refuses existing destinations, and compares Windows/VM SHA-256 values. These are transfer checks, not repository-authenticity checks. Normal-browser success is not claimed; if it also receives 403, obtain refreshed access or the five named processed files from Dryad/authors. No author message has been sent.

## 3. What has and has not been established from recordings

**Actual new recordings inspected: 0.** No measured neuron overlap, usable-animal count, copper-event count, signal missingness or baseline performance is available. Those quantities are unknown, not zero-valued biological findings. The public README schema and example are metadata evidence only. This report does not turn its example neuron count into the cohort's recorded coverage.

The offline audit is prepared to inventory normalization arrays, timestamp discontinuities, finite/constant traces, identity joins, session boundaries and descriptive reversal tables. Its per-animal, per-neuron and individual/class-by-animal outputs remain under ignored data/. Current empty table headers are placeholders, **not observed coverage**. It requires the actual author/JLD2 identity reader and processing provenance; no flattened JLD2 or nearest-position join is substituted.

The word 'original' on an exported calcium array does not certify that values were never reconstructed. Observation/reconstruction masks, timestamp units, behavioural clock mapping, radius/distance units and idx_splits semantics must be verified. Equal vector lengths cannot prove alignment. Detector saturation/artefacts cannot be diagnosed solely from numerical maxima without acquisition metadata. Left/right identities remain distinct; class-level assignments cannot silently become individual-neuron states. Calcium is not membrane voltage, and correlations are not direct anatomical coupling.

The actual author loader could not be retrieved, and no Julia executable was found. No package installations or reconstruction were started. [REPRODUCTION.md](REPRODUCTION.md) describes the exact resumption point and required private export fields. No unchecked custom JLD2 interpretation has been used to claim identity coverage.

## 4. New A–F gates and diagnostic design

See [acceptance_ledger.csv](acceptance_ledger.csv). The new requirement is **>=15 usable independent animals and the same >=15 reliably identified neuron classes in a common group of >=15 animals**. The previous eight-animal/twelve-neuron acquisition proposal is not substituted. Bilateral individual coverage must also be reported. If the broad gate fails, a narrower circuit is considered only with actual evidence and an explicitly changed scope, not an automatic whole-dataset rejection.

A, B, C, D and F remain **NOT EVALUABLE** with actual measurements. E's prospective protocol is prepared, but actual independent evaluation is **NOT EVALUABLE**. No MSE, confidence interval or biological PASS is reported. Storage/source integrity and software tests are distinct technical checks, not biological validation.

The public 27 recording identifiers were used to freeze a **15 train / 6 validation / 6 final-test** allocation before any physiological values were inspected. Exact assignments remain private; SPLIT_SUMMARY.json records counts and hash. The fixed hash-based rule is reproducible and preserves every animal's windows/sessions. Excluded animals are not moved between groups to improve performance.

[VALIDATION_DESIGN.md](VALIDATION_DESIGN.md) freezes training-only transforms, a predefined small side-resolved circuit, nested animal-held-out diagnostics within the training 15, and receiver-only AR, regularised multivariate, sensory/behaviour and low-rank/graph-constrained baselines. Validation/final-test animals do not select neurons, scaling or model parameters. Predictions use causal histories, not future smoothing. Primary behaviour is independently calibrated future velocity at boundary encounters, with reversal initiation secondary; both need real clock/unit verification and event variation. Lag/mask and acquisition-day sensitivities are prespecified. No baseline was fitted because its inputs are absent.

## 5. Mathematical design and scientific limits

[MATHEMATICAL_SPECIFICATION.md](MATHEMATICAL_SPECIFICATION.md) defines neuron/history stalks, affine multiport dynamical factors, projection restriction maps on an incidence graph, compatibility operator, masks/latent uncertainty, noise-normalized regional innovations, anchored joint-section discrepancy, predictive group coupling, a factor-based H0 filtration and circuit pushforwards. It is a proposed specification, not implemented biological evidence.

Stable Euclidean coordinates have not been verified. Named-neuron anatomy can potentially support reference graph distance or circuit membership; neither is measured micrometre separation in these moving animals. Reference electrical adjacency is not measured conductance, and unknown/partial gap annotations do not supply trustworthy negative edges.

A linear dynamical sheaf energy can be exactly a conventional state-space innovation loss. Its zero set is often consistent by construction. Added interpretability, boundary-gluing information or behavioural value must be demonstrated against conventional regional innovations and network summaries; changing the name of the loss adds no information. Unobserved neurons, common inputs, calcium filtering and behavioural feedback make direct couplings non-unique. Independent junction/circuit perturbation and rescue would be needed for direct electrical mechanism claims.

## 6. Decision and next action

**Do not implement the full sheaf yet.** There is insufficient actual-data evidence for GO, but no scientific basis for NO-GO. The remaining obstacle is precise: three processed HDF5s, the shared identity dictionary and author loader must be transferred, then audited. Storage is sufficient and the code/protocol are checkpointed. No further network retries or unrestricted dataset search are recommended.

After transfer: verify all five hashes, inspect the author loader/JLD2 structure, establish genuine observation and behavioural alignment provenance, and audit the three recordings. If they pass, obtain the remaining processed files and evaluate the >=15-animal/>=15-class gates. Only then run the frozen training-only small diagnostic. Failure of those measured gates can support a narrower PARTIAL GO or an evidence-backed NO-GO; access failure alone cannot.

Fifteen software tests pass, including exact deletion-path safeguards, symlink/hard-link refusal, invalid-download rejection, timestamp/session masks, bilateral/class distinction, genuinely shared class coverage and disjoint deterministic splits. No test is presented as a physiological validation result.
