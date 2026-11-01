# Original-recording reconciliation: decision before cohort expansion

1 November 2026. **The NWB audit does not show that all target neurons were unrecorded.** Every one of the 29 high-confidence manually assigned events has an explicit target ROI, green trace and red trace. The observed problem is incomplete fluorescence sampling, combined with a complete-case QC rule that is stricter than the authors' run-length criterion. Original export access remains blocked, so original inclusion decisions and the physical origin of missing samples remain unresolved. The nine-animal expansion is **on hold**; no further recordings were downloaded.

| Proposed experiment | Decision now | Reason / next requirement |
|---|---|---|
| Complete measured neural-state sheaf | **NO-GO under the frozen complete-case gate** | Zero complete local-state trials remain zero. Partial observations exist; recover original masks/extraction provenance before deciding whether a prospectively specified partially observed model is justified. This is not a proof that the original experiment lacked all upstream activity. |
| Experimentally constrained input-output sheaf | **PARTIAL GO for formulation and reconciliation; no fit yet** | The previous 25 manually targeted, recipient-qualified events remain a possible intervention-to-recipient design. Treat input as intended stimulation and do not infer direct local synaptic weights. Restore curation/provenance and obtain independent-animal coverage before validation. |
| WT versus unc-31 | **NO-GO pending original-data access and protocol reconciliation** | No mutant payload audited. Separate WT/unc-31 exports, independent animal matching, curation and pulse comparability are still required. |

## 1. What was recovered, and what was not

One fresh attempt at each requested OSF download failed DNS in the Codex execution sandbox. The web reader could not retrieve either link. These are execution-context failures, not evidence that the public files are absent or that normal VM/Windows access must fail. There were no repeated endpoint retries. [OSF_MANUAL_TRANSFER.md](OSF_MANUAL_TRANSFER.md) gives concise Windows download, reported version-1 checksum and SCP instructions for:

- WT: [`exported_data.tar.gz`](https://osf.io/9mecf/download), reported MD5 `7ae518fea414f664e2493238afe209bb`.
- unc-31: [`exported_data_unc31.tar.gz`](https://osf.io/ty6sp/download), reported MD5 `e5c83242baa131bca95ee58bc7b58bec`.

Names/checksums come from [Encore's reproduction record](https://github.com/catalystneuro/encore/blob/main/attempts/001075-randi-2023-01.yaml), not downloads verified here. Sizes and present HTTP availability are unverified; inspect them before transfer. The distinct `exported_data_full.tar.gz` library-compatible export is referenced by the authors' export instructions, but its direct file GUID/size have not been established.

We retrieved and hash-verified historical processing source: `pumpprobe` commit `1dbc5e0a2b609d54bc9b1c90c73d4e3bf183d3c7`, and `wormdatamodel` commit `2ab956199e3931de41a190d2b9985e961df3810c`. `source_hashes.csv`, `SOURCE_FIELDS.csv` and `fetch_sources.py` make this inspection reproducible. We did not execute the authors' fitting pipeline, deserialize pickle files, reconstruct images or fit a response classifier.

The three NWBs were re-read offline and again match their frozen SHA-256 checksums. The frozen QC directory and original NWBs were hashed before this investigation and checked again afterward. No frozen QC result was rewritten.

## 2. Reconciliation of all 29 target traces

[aggregate_failure_decomposition.csv](aggregate_failure_decomposition.csv) contains the across-event counts. Reasons below overlap and must not be added as mutually exclusive categories.

| Animal | High-confidence events | Target ROI/green/red traces found | Missing timestamps | Fails 95%-finite rule | Missing run >2 s | Low or undefined baseline SNR | Nonpositive or absent F0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 36 | 2 | 2 | 0 | 2 | 2 | 1 | 0 |
| 72 | 12 | 12 | 0 | 12 | 7 | 7 | 0 |
| 80 | 15 | 15 | 0 | 15 | 12 | 3 | 3 |
| Total | **29** | **29** | **0** | **29** | **21** | **11** | **3** |

For every event, the local event-level ledger records the manual PumpProbe ID, uniquely joined NeuroPAL label/confidence, actual signal column, finite/NaN counts in both frozen windows, longest NaN runs, baseline/noise values and all frozen reasons. The analysis verifies explicit `table_region` object references for green and red, their column-order agreement, the ROI-id/NeuroPAL joins, shared timestamps and unique target labels. No nearest-position assignment was used. The target ROI is not absent or shifted into a different column under the NWB's declared schema. This internal consistency does not independently prove the physical identity of the illuminated neuron; the original curation crosswalk remains necessary.

Both base green and base red have **identical NaN masks across their entire arrays in all three files**. The target's missing samples are not missing time coordinates, a missing entire target column, or a green-only intensity anomaly. The NWB description explicitly identifies these as signal values not inferred from imaging. Shared missingness is consistent with ROI identification/extraction/tracking loss, but we cannot distinguish original loss from conversion/masking errors without original masks, extraction logs and tracking provenance. No definitive tracking-failure or imaging-artifact cause is invented for an individual event.

A stored `InterpolatedGreenSignal` exists for each target; all 29 frozen target windows are fully finite there. Its observed-sample values agree exactly with base green. It is the **same ROI under a processed representation**, not a newly found measured trace. We neither generated interpolated values nor used the stored fills to pass the primary gate. Red fluorescence is a reference measurement, not an alternative calcium response. There is no stored matchless trace series in these NWBs. Their analysis groups are empty; there are no saved `resp_neurons_by_stim`, authors' exclusion codes or classifier/autoresponse flags to consult.

None of these 29 target windows contains a negative sample. Subject 36's previously documented extreme negatives remain real QC concerns elsewhere, but they do not explain this target-trace failure set. Its two events fail completeness/long-gap rules; one also lacks enough adjacent baseline samples to estimate noise. No values were clipped or offset.

## 3. Completeness: a controlled comparison, not a changed gate

Our frozen primary rule requires at least 95% finite values in **each** -10–0 s and +1–20 s window, plus a missing-run limit, baseline and other QC checks. The authors' `Fconn.nan_ok` instead checks whether the **longest contiguous run** of missing values exceeds a fraction of an analysis interval. Scattered missing samples can fail the former while passing the latter. The pipeline normally interpolates fluorescence but retains an original `Signal.nan_mask` for this check. These criteria are not interchangeable. Sources: [Fconn.nan_ok](https://github.com/leiferlab/pumpprobe/blob/1dbc5e0a2b609d54bc9b1c90c73d4e3bf183d3c7/pumpprobe/Fconn.py#L1751), [Signal preprocessing](https://github.com/leiferlab/wormdatamodel/blob/2ab956199e3931de41a190d2b9985e961df3810c/wormdatamodel/signal/signal.py#L104).

With the **same frozen windows**, a five-percent longest-run-only component passes 0/3/3 events in animals 36/72/80: six total. Five events fail **only** the original total-completeness rule. These component diagnostics do not replace that rule or certify a complete local-state trial.

Additional timestamp-window diagnostics reproduce only the run-length calculation on NWB missingness, not the full original workflow:

| Animal | Five-percent run limit, -30 to +30 s | Five-percent run limit, -30 to +60 s | Thirty-percent run limit, -30 to +60 s |
|---|---:|---:|---:|
| 36 | 1/2 | 1/2 | 2/2 |
| 72 | 5/12 | 6/12 | 10/12 |
| 80 | 6/15 | 5/15 | 10/15 |
| Total | **12/29** | **12/29** | **22/29** |

Original stimulus-volume indices, actual `Dt`, source masks and matchless substitutions are unavailable. These are timestamp-based component checks on clipped NWB windows where needed; they are **not author inclusion decisions**. No favourable window was chosen as a replacement.

A reproducible source discrepancy also matters: the pinned pipeline requests `--nan-th:0.05`, but its `responses_detect.py` parser does not handle that option and leaves `nan_thresh=0.3`. Unit tests confirm this from the AST. Thus documented five-percent intent, class defaults, example CLI and effective pinned script settings must not be collapsed into one claimed original configuration. Saved per-recording settings or an author clarification are required. Source: [pipeline](https://github.com/leiferlab/pumpprobe/blob/1dbc5e0a2b609d54bc9b1c90c73d4e3bf183d3c7/pipelines/fit_kernels_multi_datasets.sh), [parser](https://github.com/leiferlab/pumpprobe/blob/1dbc5e0a2b609d54bc9b1c90c73d4e3bf183d3c7/scripts/fconnectivity/responses_detect.py#L23).

The earlier zero count was correct **for our stated complete-case definition**, but cannot establish that the authors lacked measurable target activity or would have rejected all 29 events. That broader interpretation would be too strong. Conversely, dense interpolated export arrays cannot by themselves establish recovered measurements. A prospective partially observed formulation may be defensible, but cannot be labelled a pass of the frozen complete-state gate.

## 4. Autoresponses and original exclusions: exact information boundary

Nineteen events have at least five actually observed samples in both frozen pre/post windows and a positive baseline, allowing an observed-sample median change to be described: 2/8/9 events in the three animals. These descriptive values are retained locally. Sparse fluorescence changes are **not** a demonstrated autoresponse: they do not implement the original amplitude/derivative-duration classifier, preprocessing or curation, and missingness can bias them.

The [published inclusion procedure](https://doi.org/10.1038/s41586-023-06683-4) normally requires a stimulated-cell autoresponse; it describes preprocessing and a special RID analysis exception. None of the 29 events here targets RID, so that exception does not certify any of them. The authors' code stores detected response ROI lists in `resp_neurons_by_stim`; target membership defines its autoresponse flag. `targeted_neuron_hit` is a separate flash/bleaching crosscheck. `stim_neurons` can combine tracking-based closest-neuron assignment, manual overrides and later exclusions. Positive export indices alone do not prove a manual assignment or full inclusion. Missing values may arise when a neuron is not identified in a volume; matchless image extraction is a distinct alternative to temporal interpolation. Sources: [Fconn response/target methods](https://github.com/leiferlab/pumpprobe/blob/1dbc5e0a2b609d54bc9b1c90c73d4e3bf183d3c7/pumpprobe/Fconn.py#L655), [manual overrides](https://github.com/leiferlab/pumpprobe/blob/1dbc5e0a2b609d54bc9b1c90c73d4e3bf183d3c7/pumpprobe/Fconn.py#L546).

Accordingly, **original retained/excluded status, original detected autoresponse and full published inclusion are UNKNOWN / NOT EVALUABLE for all 29 events here**. Unknown is not false, and it is not a new exclusion. A classifier run on imputed NWB data would not recover the original decisions.

The flat export writer provides `*_gcamp.txt`, `*_t.txt`, `*_labels.txt`, `*_stim_volume_i.txt`, `*_stim_neurons.txt` and `*_ds_name.txt`. It does **not** write original NaN masks, matchless provenance or response-list flags. Target code -2 indicates an export exclusion but can arise for different reasons; -1 and -3 need separate interpretation. The full library export or safe author-exported state is needed for complete decision reproduction. Source: [Funatlas.export_to_txt](https://github.com/leiferlab/pumpprobe/blob/1dbc5e0a2b609d54bc9b1c90c73d4e3bf183d3c7/pumpprobe/Funatlas.py#L279).

## 5. Stable original-export matching and resumable workflow

`inspect_exports.py` implements bounded, plain-text tar inspection with no extraction or pickle execution. Matching uses acquisition keys **20220331_114306**, **20230309_143050**, **20230314_164209** from source dataset names, not the integer OSF export ordinal. No actual cross-repository match is claimed until those files arrive.

The scanner requires unique session matches, compatible frame/event counts, stimulus-volume alignment and unique label joins. Disagreements are recorded rather than repaired by signal correlations or arbitrary index shifts. It verifies the reported version-1 archive checksum before accepting a known filename, additionally hashes SHA-256, and records the inventory locally. `_t.txt` is generated from zero and a fixed `Dt`, whereas NWB times have a nonzero origin: source indices and time origins must be reconciled explicitly. Dense export traces remain measured-versus-imputed **unresolved** without masks. The code has a synthetic archive integration test; it has not been exercised on actual OSF payloads.

The local event-level CSV is:
`data/worm_data_access_2026/original_reconciliation_v1/event_failure_decomposition.csv`.
It contains all 29 event rows, explicit unknown original statuses and exact NWB failure reasons. It is intentionally **not published to GitHub**. Only across-event aggregates, code, source provenance and instructions are published.

## 6. Two model formulations and fair baseline comparisons

**A. Measured local neural dynamics.** Restrictions would relate reliably observed upstream/downstream neuronal states and their delays on anatomically constrained edges. It requires appropriate observed-history coverage, credible neuron identity and a calcium observation model. The current complete-state gate fails; missingness is not solved by giving latent reconstructions the status of measurements. A prospectively specified model with explicit observation masks is a different formulation and should be tested for practical identifiability before mechanistic claims. Anatomy does not make local weights uniquely identifiable from effective calcium propagation.

The appropriate comparisons are receiver-only autoregressive/history forecasting, an intrinsic baseline with direct stimulus effects, and a connectome-constrained dynamical baseline such as Creamer–Leifer–Pillow, all using the **same causal history and observed evaluation points**. Future-data smoothing must not be compared with causal forecasting. Anatomically matched shuffles and observation-mask-preserving controls test the contribution of neighbours. Current data cannot supply a replicated direct-edge animal-held-out test; no baseline was fitted or evaluated numerically.

**B. Intervention-to-population response.** Restrictions would relate an explicitly recorded intended stimulus, its target identity/history and receiver baseline to an observed population-response window. Such effective, time-dependent maps can include multiple synaptic steps and extrasynaptic propagation; they are not direct local synaptic restrictions. The present 25 recipient-qualified events remain a design foothold. If target activation is unverified, the estimand is response to **intended stimulation**, not confirmed stimulated-cell activity. Choosing an intention-to-stimulate analysis would be an explicit new estimand, not a disguised relaxation of author inclusion rules.

Compare against a strong receiver-history predictor, an identity-blind stimulus-timing response model and a training-only target-specific empirical response baseline. Match stimulus times, baseline states, receiver identity, observation masks and target frequencies in nulls; group all versions of an animal together. Fit transforms/maps only on training animals. Require the prespecified held-out improvement and animal-level uncertainty, and distinguish original autoresponse-selected analyses from a prospective intention-to-stimulate analysis. No model or compatibility metric was fitted here, so these are required future comparisons, **not successful validation results**.

## 7. Expansion decision and exact author enquiry

**Do not download the nine additional NWBs yet.** The expansion condition—recovery of reliably measured target traces and trustworthy event identities from originals—has not been established. The current result is not that nine more animals would be useless forever; it is that increasing the same representation cannot resolve an untested mask/conversion/curation issue.

After WT transfer, first reconcile these three animals. If exports or matchless products recover measurement-backed traces with clear masks and event identities, preserve the original zero-complete-case result, document recovered representations separately and then execute the existing deterministic nine-animal QC manifest. Keep all selected animals, including failures, in the ledger. If the originals are only dense interpolated arrays, or originals confirm genuinely unmeasured upstream intervals, do not call that recovery or proceed merely to inflate n. Prioritize an explicitly partially observed or input-output formulation, or obtain better measurements.

If flat exports do not settle the issue, request a **small processed package for these three named sessions**, without raw video:

1. Original green/red arrays and original `nan_mask`, alternate matchless arrays, and provenance distinguishing image extraction from temporal fills.
2. Source `Dt`, volume timestamps, stimulus-volume indices, tracking-to-reference IDs, NeuroPAL labels/confidences and manual target overrides/comments.
3. Safe JSON/CSV of `stim_neurons`, `resp_neurons_by_stim`, `targeted_neuron_hit`, `i0s`, `i1s`, `shift_vols`, bubble/manual response flags and the actual per-session detector settings.
4. Clarification of the ignored `--nan-th` option versus actual saved settings, and original pulse metadata for WT/unc-31 comparisons.

These files could resolve conversion versus extraction loss, curation and original classifier status. No author message was sent. An author enquiry is not yet necessary for downloading advertised open files, but may be necessary for original masks, matchless provenance and decision state omitted from the text exports.

## 8. Reproduction and preservation

Run `bash review/worm_data_access_2026/original_reconciliation_v1/run_all.sh` from the repository root. It uses the existing Python 3.12/h5py environment without altering it, verifies cached source hashes and executes 12 new tests. Optional `--wt-archive` invocation is in the transfer guide. On a fresh checkout, obtain the three original NWBs and use the supplied metadata bootstrap; no observations are part of this commit.

Tests cover contiguous-run versus total-completeness distinctions, exact boundary behaviour, no interpolation, stable-session matching, ambiguous sessions, bounded archive reads, unsafe paths, version mismatch, a synthetic export crosswalk and the source CLI discrepancy. All 12 pass. Frozen QC code/rules and prior tests remain intact. A source-only dependency allowlist is included in the GitHub publication so this reconciliation can run without publishing the prior event/ROI checkpoints. No raw or individual-event observations, source cache, NWBs, archives, fitted model or higher-level metric is pushed. The final decision is constrained by observed missingness and explicitly missing originals, not by an assumption that all target activity was lost.
