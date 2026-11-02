"""Compile aggregate decision from completed audits; never fit a model."""
import csv, json, hashlib
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parent
P=R.parents[1]
def main():
 d=pd.read_csv(R/'liberal_detection_radius_QC_sensitivity.csv')
 s=d.groupby('programme').agg(max_receptor_fraction=('competent_fraction','max'),max_local_candidate_fraction=('candidate_fraction','max'),max_local_candidate_count_in_one_tissue=('spatial_candidate_receivers','max')).reset_index()
 s.to_csv(R/'liberal_coverage_bounds_summary.csv',index=False)
 sub=pd.read_csv(R/'all_fibroblast_subtype_receptor_bounds.csv')
 bounds=sub.groupby('programme').necessary_upper_bound_gate.sum().to_dict()
 nom=pd.read_csv(R/'diagnostic_receiver_coverage.csv')
 nom=nom[(nom.population=='author_fibroblast_labels_conf_ge_08')&(nom.radius_multiplier==2.5)]
 nom.to_csv(R/'nominal_coverage_by_specimen.csv',index=False)
 totals=nom.groupby('programme')[['eligible_cells','raw_receptor_ge1','screened_receptor_competent','immune_senders','exposed']].sum()
 totals.to_csv(R/'nominal_coverage_summary.csv')
 frozen=json.loads((R/'BIOLOGICAL_FROZEN_INPUTS.json').read_text())
 for name,digest in frozen.items(): assert hashlib.sha256((P/name).read_bytes()).hexdigest()==digest,name
 rows=[
 ('A1','PASS','100% unique expression-row and composite specimen/cell joins; all 2,131,098 centroids finite. Tissue boundaries respected; no exact tissue outline.'),
 ('A2','PASS','Unscaled Xenium x/y centroids and cell areas; standard units micrometres and square micrometres; coordinates match obsm/spatial.'),
 ('A3','NOT EVALUABLE','Broad labels assigned to 100%; credible independent validation unavailable. Refined QC label completeness 78.8195%; annotation agreement is not accuracy.'),
 ('A4','NOT EVALUABLE','Only PDGFRA of six frozen structural markers measured; >=2-marker eligibility cannot be assessed. FAP/PDPN diagnostics cannot replace gate.'),
 ('B1','NOT EVALUABLE','IFNG, TNF, IL1B ligands and canonical receptor genes present; external perturbation payload and gene-disjoint dictionary not inspected.'),
 ('B2','NOT EVALUABLE','No verified cytokine signatures: >=10 response genes and >=5 detected in >=20% receivers cannot be evaluated.'),
 ('B3','NOT EVALUABLE','Culture replicate effects, significance and reproducibility unavailable; four cell-line cultures would not be independent tissue patients.'),
 ('C1','FAIL','Across all tested author fibroblast subtypes, IFNG and IL1B each have zero QC/subtype scenarios satisfying even receptor-only >=100 and >=20% in >=3 tissues; two programmes cannot pass.'),
 ('C1 patient identity','NOT EVALUABLE','Workbook supports 16 tissues/8 provisional pairs; authoritative S1 unavailable. No patient-grouped statistics performed.'),
 ('C2','NOT EVALUABLE','Aggregate negative-control screen evaluated, but per-probe calibrated controls and protein occupancy absent. Liberal bounds already fail.'),
 ('C3','PASS','Coverage sensitivity across 1.5/2.5/4 diameters, raw thresholds 1/2/3 and QC count floors 20/50/100/200; all canonical subunits retained. This methodological PASS does not pass C.'),
 ('D1-D4','NOT EVALUABLE','No response prediction fitted: prerequisites failed/unverified. No MSE improvements or patient confidence interval claimed.'),
 ('E1-E4','NOT EVALUABLE','Coverage sensitivity completed; predictive sign stability and downsampling not evaluated because no model fitted.'),
 ('F','NOT EVALUABLE','Eight workbook patients provisionally recovered, not seven slide categories; S1 not verified. Affected/adjacent pairing is not independent case/control replication.')]
 with (R/'biological_criterion_assessment.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['criterion','status','evidence']);w.writerows(rows)
 source=P/'data/crohn_xenium/scp3259/xenium5k_scVI_CellCharter_batch2_reannot_SCP.h5ad'
 inv={'file':source.name,'bytes':source.stat().st_size,'sha256':'74759d2d474964dcb537533127f0683bdb9180030f3b94a8f3f8271a061494bb','cells':2131098,'genes':5001,'format':'H5AD / CSR raw counts and log-normalised X','source':'https://singlecell.broadinstitute.org/single_cell/study/SCP3259','private_linkage_published':False}
 (R/'SOURCE_INVENTORY.json').write_text(json.dumps(inv,indent=2)+'\n')
 private=P/'data/crohn_xenium/biological_validation_v1/IMMUTABLE_SOURCE_MANIFEST_PRIVATE.json'
 if not private.exists():
  manifest=dict(inv,workbook=json.loads((P/'data/crohn_xenium/biological_validation_v1/workbook_linkage_PRIVATE.json').read_text()))
  private.write_text(json.dumps(manifest,indent=2)+'\n');private.chmod(0o444)
 decision={'decision':'NO-GO','scope':'Validated local signalling and multiscale sheaf with current measured assay under frozen criteria','spatial_ingestion':'PARTIAL GO','prediction_fitted':False,'patient_mapping_S1_verified':False,'external_signatures_verified':False,'subtype_receptor_bound_passing_scenarios':{k:int(v) for k,v in bounds.items()}}
 (R/'DECISION.json').write_text(json.dumps(decision,indent=2)+'\n')
 lines=['| Programme | Maximum complete-receptor fraction | Maximum local-candidate fraction |','|---|---:|---:|']
 for x in s.itertuples(): lines.append(f'| {x.programme} | {100*x.max_receptor_fraction:.3f}% | {100*x.max_local_candidate_fraction:.3f}% |')
 table='\n'.join(lines)
 report=f'''# Crohn Xenium biological validation — final decision, 2026-11-02

**NO-GO for validated local signalling or multiscale sheaf construction with this measured assay under the frozen acceptance criteria.** PARTIAL GO for descriptive atlas use and provenance completion. This conclusion uses the actual downloaded spatial measurements, not publication metadata. It does not establish biological absence of signalling. No model, disease classifier, compatibility energy or persistent homology was fitted/computed.

## Actual acquisition and integrity

The public SCP3259 H5AD is {inv['bytes']:,} bytes; SHA-256 `{inv['sha256']}`. It contains **2,131,098 cells × 5,001 genes**, 345,626,891 sparse raw-count entries. Every stored count was checked: zero negative, noninteger or nonfinite values. Original files remain unchanged. See SOURCE_INVENTORY.json and the private read-only source manifest. The earlier transfer checkpoint is preserved under stage1_checkpoint/.

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

{table}

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
'''
 (R/'REPORT.md').write_text(report)
 print(s.to_string(index=False));print(totals.to_string())
if __name__=='__main__':main()
