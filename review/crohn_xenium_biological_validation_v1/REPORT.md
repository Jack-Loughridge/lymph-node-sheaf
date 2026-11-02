# Crohn Xenium biological validation — final decision, 2026-11-02

**NO-GO for validated local signalling or multiscale sheaf construction with this measured assay under the frozen acceptance criteria.** PARTIAL GO for descriptive atlas use and provenance completion. This conclusion uses the actual downloaded spatial measurements, not publication metadata. It does not establish biological absence of signalling. No model, disease classifier, compatibility energy or persistent homology was fitted/computed.

## Actual acquisition and integrity

The public SCP3259 H5AD is 8,984,620,838 bytes; SHA-256 `74759d2d474964dcb537533127f0683bdb9180030f3b94a8f3f8271a061494bb`. It contains **2,131,098 cells × 5,001 genes**, 345,626,891 sparse raw-count entries. Every stored count was checked: zero negative, noninteger or nonfinite values. Original files remain unchanged. See SOURCE_INVENTORY.json and the private read-only source manifest. The earlier transfer checkpoint is preserved under stage1_checkpoint/.

All expression-row IDs and composite (specimen, barcode) joins are unique. Individual barcodes alone repeat 527 times: do not join them globally. All centroids are finite; obsm/spatial equals source x/y. Standard [Xenium output units](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/1.9/tutorials/outputs/xoa-output-understanding-outputs) are micrometres and cell areas square micrometres. There are 16 tissue categories, no z coordinates or exact tissue outlines. CellCharter neighbourhoods are computational clusters, not independently validated anatomical regions.

Both raw CSR counts and normalised CSR X are available. In 16,000 sampled rows, after aligning independently ordered sparse indices, X matches log1p(10,000 × count / source total_counts). Exported panel sums differ from source total_counts for 1,458,553 cells; the difference is not explained by exported negative-control counts. The upstream count universe is unresolved. QC therefore uses actual exported count sums, never assumes source totals equal panel totals.

Dense float32 atlas storage alone would require 42.6 GB. Initial free disk was 9.72 GB and available RAM about 13.2 GiB. Sparse 25,000-row streaming and a 38-gene cache avoided dense conversion, duplicate downloads and alteration of the existing analysis environment.

## Samples and identities

The private workbook was inspected and hashed. Full tissue identifiers, embedded biobank identifiers where needed and workbook rows provisionally resolve **16 tissues into eight paired patients**. The seven source donor_id categories are slide identifiers. Patients 5 and 6 share a slide and are not collapsed. Exact linkage remains under ignored data/; public specimen aliases disclose no patient/biobank mapping.

**Supplementary Table S1 remains unavailable here, so the mapping is provisional and no patient-grouped analysis was performed.** No ambiguity was silently resolved using slide ID alone. Populated inflammation_status values agree with explicit specimen-name status, but 1,484,977 cells have blank source status. Reports distinguish name-derived status; adjacent is not described as healthy. All eight provisional pairs are represented; authoritative confirmation is pending.

Primary QC requires >=100 exported counts, >=20 detected genes, valid coordinates and positive area. **1,425,592 cells (66.895%) pass**, with tissue fractions 36.38–83.48%. All receive broad source labels. Refined annot_new labels are missing for 667,030 cells, including 301,948 QC cells: **78.8195% QC refined-label coverage**. Among annotated QC cells, mapped broad/refined agreement is 89.5015%; this is concordance between related annotations, not independent accuracy. Missing categorical code -1 is explicitly handled and tested.

Only PDGFRA of the six frozen structural fibroblast markers is measured. The >=2-marker eligibility screen is **NOT EVALUABLE**, not evidence of absent fibroblasts. Post-panel FAP/PDPN diagnostics do not replace this screen. Consensus author-fibroblast/broad-mesenchymal labels identify 101,892 QC cells, including 58,350 with source confidence >=0.8; 4,754 carry the inflammatory fibroblast label. Confidence is not independent identity corroboration. FAP or PDPN is detected in 27,387 consensus-labelled fibroblasts, both in 3,821. Broad mesenchymal cells total 408,091 QC cells. Independent marker/protein/morphology validation remains incomplete.

## Measured signalling coverage

All IFNG/IFNGR1/IFNGR2, TNF/TNFRSF1A (also TNFRSF1B), and IL1B/IL1R1/IL1RAP genes are measured. Complete complexes require co-detection in the same cell. OSMR is absent: canonical primary OSM analysis is excluded. LIFR/IL6ST inventory is not a validated alternative response map or a rescue of the primary analysis.

The fixed aggregate-control screen requires each ligand and receptor subunit to exceed per-cell negative controls and reach max(2, 1 + the tissue QC control-count 99.9th percentile). Floors are two counts in all tissues. This is a conservative aggregate-control screen, not calibrated per-probe background inference. Senders are immune cells; searches stay inside each tissue. Radii are 1.5, 2.5 and 4 median segmentation-equivalent diameters (diameters 6.77–8.85 micrometres). Truncated tissue-edge exposure cannot be fully assessed without outlines.

At nominal 2.5 diameters in the 58,350 confident consensus-labelled fibroblasts, screened complete receptors occur in IFNG 12, TNF 1,031 and IL1B 7 cells. Spatial immune-ligand candidates are **IFNG 0, TNF 4 and IL1B 0**. These are transcript-neighbour candidates, not measured cytokine binding or causal effects.

A permissive sensitivity removes background correction, requires only >=1 molecule for all receptor subunits and ligands, and evaluates QC floors 20/50/100/200, all three radii and four author-label populations. The maxima below are across tested tissues/settings, not a selected passing model. Counts >=2 and >=3 are also retained in the full sensitivity table.

| Programme | Maximum complete-receptor fraction | Maximum local-candidate fraction |
|---|---:|---:|
| IFNG | 3.196% | 0.935% |
| IL1B | 1.539% | 0.344% |
| TNF | 26.876% | 3.316% |

More decisively, an exhaustive audit of every author fibroblast subtype at QC floors 20/100/200 and confidence >=0.8 **ignores ligand, background and distance altogether**. IFNG and IL1B each have **zero subtype/QC scenarios** with >=100 complete-receptor cells and >=20% receptor competence in >=3 tissues; TNF has 11 scenarios. This is a necessary upper bound, not exposure success. It prevents a broad receiver denominator alone from explaining failure. Resolving patient pairing cannot create missing measured receptor subunits; definitive patient-level statistics remain NOT EVALUABLE.

Accordingly, two primary programmes cannot meet the exposure gate in the audited measured stromal populations. This is an assay-resolution limitation: dropout, persistent receptor protein, secretion and diffusion may sustain real signalling despite sparse RNA. No receptors were imputed, subunits dropped or spatial radii inflated to obtain a pass.

## External perturbations and acceptance ledger

Project searches found no GSE307526 processed cytokine payload, response DEG supplements or authoritative S1. Bounded NCBI/publisher requests failed DNS/access in this executor. This does not prove public access or the user's VM/browser is blocked. [EXTERNAL_INPUTS.md](EXTERNAL_INPUTS.md) lists exact small filenames, primary links and Windows/SCP transfer commands. No repeat full-atlas download is needed.

The published experiment uses a stromal cell line, not matched spatial patients or an organoid donor replication design. IFNG, TNF and IL1B external effect sizes, signed significant targets, culture reproducibility, >=10-gene panel overlap and actual >=5-gene receiver detection remain **NOT EVALUABLE** until those payloads are inspected. Cell-state DEGs in H5AD uns are not cytokine signatures. Four culture replicates cannot be counted as independent patients.

See biological_criterion_assessment.csv for every gate: A1/A2 pass technical checks; credible independent annotation A3/A4 and response gates B remain not evaluable; C1 fails measured coverage; predictive D/E and independent-patient F remain not evaluable. A methodological C3 sensitivity pass is not a biological C pass. No held-out MSE, confidence interval or predictive success is claimed. COMMOT was not installed or fitted after its prerequisites failed; KD-tree coverage enumeration is not COMMOT.

## Next experiment and stopping decision

Do not expand signalling modelling or construct the full multiscale sheaf from these measurements. Retrieve the five small provenance/perturbation files to close mapping and response-dictionary gaps; they will not by themselves fix observed receptor coverage. Keep the atlas for descriptive questions only.

A defensible next local-validation experiment requires deeper relevant transcript coverage or independently measured receptor/cytokine protein and physiological responses, such as stromal pSTAT1 and NF-kB activity, with matched perturbations and negative controls. Validate transfer from cultured cells explicitly. Require the same frozen coverage and independent-response gates before modelling.

If an improved assay passes, freeze a signed external response dictionary excluding ligand, receptor and annotation targets; retain whole affected/adjacent patient pairs in outer folds. Fit normalisation, calibration and radius choices using training patients only. Compare against a target-excluded receiver-intrinsic model, matched within-tissue spatial/cell-type sender nulls and simpler neighbourhood models. Require >=5% held-out MSE improvement versus both primary baselines for two programmes, patient-level uncertainty excluding zero, and prespecified radius/downsampling stability. Eight pairs warrant careful uncertainty analysis and external replication, not automatic confidence. No assumption that affected tissue is less compatible is warranted.

## Reproduction and confidentiality

Run bash review/crohn_xenium_biological_validation_v1/run_biological_audit.sh from the project root, then run finalise_report.py using the audit Python. Tests cover complete complexes, missing-category codes, sparse offsets, spatial bounds and ambiguous/shared-slide joins. Original mathematical specification and biological thresholds are hash-verified unchanged. All exact linkage, raw matrices, cell caches and any future cell predictions stay under ignored data/. Only code and aggregate reports are eligible for publication.
