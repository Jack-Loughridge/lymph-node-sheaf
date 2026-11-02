# Spatial-context-to-stromal-response feasibility — checkpoint

**Decision: BLOCKED / NOT EVALUABLE for the alternative representation; original Crohn NO-GO remains frozen.** A separate prospective protocol and a measured-data background audit are complete. No new response score, predictor, sheaf or higher-level metric has been fitted or calculated. This is not a successful alternative model, and blocked data access is not a scientific failure of the alternative.

## Background discrepancy

The previous pooled per-cell control screen is **not equivalent** to 10x's per-control-probe background metric. It combines probe and codeword controls, omits the number of control targets and imposes a minimum-two count rule. That can suppress real low-count observations. The official metric is a mean background estimate, not a universal rule for gene/cell positivity. See [method audit](BACKGROUND_METHOD_AUDIT.md) and [10x documentation](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/1.9/analysis/xoa-output-analysis-summary).

Actual exported counts total 6,745 negative-control probes and 5,197 negative-control codewords across 2,131,098 cells. Per-control metrics can be calculated only conditionally until the run's control-target denominator is verified. Individual control features are absent from this H5AD. The new audit changes no original exposure calls or criterion statuses; all 53 recorded original files/specifications were hash-verified unchanged. The original >=1-count complete-receptor sensitivity did not use the disputed threshold, so the discrepancy alone does not overturn that finding.

## A distinct estimand

[PREREGISTRATION.md](PREREGISTRATION.md) asks whether local measured immune-ligand context improves **patient-held-out prediction** of externally defined stromal response scores beyond intrinsic expression, cell composition, independent anatomy and section-wide ligand abundance. Receiver-type competence requires independent evidence; each receiver need not show receptor mRNA. This cannot retrospectively pass the original receptor-exposure gate.

The prospective relation is an affine neighbourhood-context/receiver-state factor, with error calibrated from training patients. Its biological interpretation is conditional association, not direct causal ligand signalling or identified pairwise coupling. Context is an aggregate over neighbouring senders, not a known direct physical edge. No sheaf has been constructed.

The protocol fixes IFNG/TNF/IL1B, external DEG selection, gene disjointness, target-excluded normalisation, primary radius and sensitivity grid, eight patient-pair outer folds, nested training-only model selection, strong nonlinear and linear baselines, composition/anatomy controls, matched spatial nulls and patient-level uncertainty. It retains >=5% improvement and >=2 programmes, plus coverage and stability gates, while explicitly defining the new context-support gate separately from the original receptor gate. Fine inflammatory subtype labels and expression-derived neighbourhoods cannot be used as independent confounder controls if they encode target genes.

## Epigenetic persistence limits interpretation

The study reports an inflammatory fibroblast state that persists after withdrawal/blockade of its inducing cytokine mixture, together with chromatin-state changes and effects of epigenetic modifiers. This supports a history-dependent explanation for some states. The mixture involved IL1B, TNF, IL17A and OSM; it does not establish that IFNG alone induces persistent states. [Original study](https://www.nature.com/articles/s41590-026-02617-0).

Consequently, a static neighbourhood may predict an anatomical/history-associated inflammatory state without identifying its current causal driver. Conversely, a competent persistent fibroblast can retain a response programme after local inducing signals disappear. Neither a positive nor a negative spatial association distinguishes these possibilities. Current data cannot reconstruct exposure history or prove causal ligand dependence. Eight held-out patients also provide limited uncertainty resolution; fold-loss bootstrap intervals are not substitutes for external replication.

## Exact remaining blockers

The local supporting-file search is negative. One bounded attempt per named source file failed DNS in this executor; publisher web retrieval also failed. The runtime plugin's environment-status tool and network-policy file were unavailable here, so no VM-wide networking conclusion is drawn. File endpoints and access errors are in supporting_file_access.json.

1. Authoritative Supplementary Table S1, to verify the private eight-patient mapping, including the shared-slide pair.
2. Three external cytokine DEG tables and GSE307526 replicate expression metadata, to verify directions, effects, target overlap and reproducibility.
3. Independent receptor-competence evidence for the eligible receiver type and two programmes; the cell-line dictionary alone cannot establish competence of every tissue subtype.
4. Independently drawn anatomy/section outlines or a valid external anatomical reference. Current CellCharter labels do not satisfy this independently.
5. Run-specific control-target denominator and negative-control feature profiles for calibrated neighbourhood background. Separate aggregate control columns suffice for a global mean with a verified denominator, not for probe heterogeneity or spatial control distributions.

[MANUAL_INPUTS.md](MANUAL_INPUTS.md) contains exact small-file links and transfer commands. No full-atlas redownload is needed. Access must be resolved before patient-held-out prediction; the current checkpoint does not claim these biological gates passed.

## Reproduction and decision rule

Use the existing h5py/scipy environment without modifying the microbiome environment:

```bash
python3 review/crohn_spatial_context_preregistered_v1/retrieve_supporting_files.py
/home/Jack/micromamba/envs/humann_py312/bin/python review/crohn_spatial_context_preregistered_v1/background_audit.py
/home/Jack/micromamba/envs/humann_py312/bin/python -m unittest discover -s review/crohn_spatial_context_preregistered_v1 -p test_protocol.py -v
/home/Jack/micromamba/envs/humann_py312/bin/python review/crohn_spatial_context_preregistered_v1/protocol_checks.py
```

Nine tests pass: original/protocol hash integrity, pair separation, S1 requirement, duplicate/ambiguous units, feature exclusion, correct control denominator and error-page rejection. The readiness check fails closed and does not launch fitting. PREREGISTRATION_SHA256.json records the protocol/code before any new outcome evaluation. This is repository preregistration, not an external registry entry. All exact linkage and recordings remain private.

If the required files become available, audit them and run only the preregistered experiment if its gates pass. If the properly controlled representation fails, stop using this snapshot for the required biological validation and move to direct physiological perturbations with measured stimulus timing, receiver dynamics, controls and independent animals/cultures. That would establish inadequacy of this assay/representation for the intended validation, not impossibility for all static spatial transcriptomics. No further arbitrary receptor-threshold adjustments are recommended.
