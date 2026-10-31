# Prespecified biological adequacy gates — spatial signalling sheaf (2026-10-31)

These are **research-stage decision thresholds**, NOT empirically established
universal standards. They must NOT be lowered after inspecting results.
Predeclare candidate pathways, thresholds, training/validation units and
parameter grids before testing. Preserve the frozen initial pilot.

## Scope of the decision

- **TECHNICAL SHEAF:** A mathematically valid factor-graph sheaf already
  exists. Passing algebraic tests alone is not a biological success.
- **BIOLOGICAL-RESOLUTION GATE (A–C):** Sufficiently resolved, spatially
  meaningful local signalling relationships for an interpretable exploratory
  analysis. Failure triggers a source/data audit or dataset pivot.
- **BIOLOGICAL-VALIDATION GATE (D–E):** Evidence that signalling information
  explains independently observed response beyond intrinsic and null models.
  Passing is required **before** treating sheaf energy or multiscale features
  as a validated biological organisation metric.
- **CLINICAL/DISEASE GATE (F):** Independent donor replication. No clinical
  interpretation of a four-specimen pilot is permitted.

## A. Cell identification and spatial measurement

1. >=95% of expression-passing cells must have unambiguous, unique metadata
   joins and valid spatial coordinates; record tissue/FOV boundaries.
2. Verify conversion to physical distance (micrometres) before claiming
   biologically comparable spatial scales.
3. At least 90% of QC cells must receive a credible **broad lineage** label
   (epithelial/immune/stromal/endothelial/other) from author-provided labels
   or external reference transfer, with independent marker/protein/morphology
   corroboration. At least 70% of *epithelial cells* must have sufficiently
   supported subtype labels for the chosen response dictionary. High model
   posterior without validation does not count.
4. Record a marker/protein and anatomy concordance audit. A label-proxy-only
   result is exploratory and does not pass A.

## B. Independent response-map resolution

1. >=2 mechanistically distinct programmes (aim 3: IFNG, TNF, BMP), with
   canonical complete receptor components and separately measured ligand
   genes, must be identifiable without input/response gene reuse.
2. Each passing programme needs >=10 experimentally responsive, non-input
   genes on the spatial panel; >=5 of these must be detected in >=20% of
   correctly typed eligible receiver cells.
3. Check source organoid donor reproducibility using the published
   matched-donor permutation criterion (>95th percentile of their null, if
   available). If it cannot be reconstructed, state this explicitly and
   seek independent perturbation evidence. Avoid adding overlapping pathways
   as though their response targets were independent.

## C. Spatial signalling coverage

1. For each of >=2 independent passing programmes, require >=100 receptor-
   competent **exposed** receivers across >=3 independent specimens, with
   >=20% of appropriate, confidently typed receivers exposed in each of
   >=3 specimens. Report both absolute counts and fractions.
2. An exposure must exceed a prespecified noise/background floor; merely
   being numerically >0 is not sufficient. Distinguish genuinely missing
   signal from technical dropout and truncated FOV exposure.
3. Do not achieve coverage by inflating a radius, dropping receptor
   subunits, or raising background signal. Compare plausible ranges and
   explicit positive/negative controls, including no-receptor/null routes.

## D. Held-out biological predictive utility

1. Freeze external signed response maps and tune calibration only on
   training specimens. Use grouped donor-/specimen-held-out prediction
   (never random cells/FOVs from the same specimen across folds).
2. Require >=5% lower held-out response MSE relative to BOTH:
   (i) a strong receiver-intrinsic, target-excluded predictor; and
   (ii) a spatially and cell-type-matched sender-shuffle/null model.
   The type-intercept-only baseline is insufficient.
3. The improvement must hold in >=2 mechanistically distinct programmes,
   and the 95% uncertainty interval on donor-level paired differences
   must exclude zero. If there are too few donors for meaningful intervals,
   mark D **not evaluable**, never as passed.
4. Report performance separately for truly exposed receivers and all
   eligible receivers; evaluate FOV-edge bias and gene-panel missingness.

## E. Stability and resolution invariance

1. Across biologically plausible spatial ranges (e.g. 1.5/2.5/4 median
   cell diameters after physical calibration), the sign of the held-out
   added-value estimate must agree in >=2/3 radii.
2. Under 20% cell downsampling with appropriate region/density matching,
   retain the effect's direction in >=80% of repetitions and at least half
   of the original effect magnitude in the median repeat.
3. Repeat at more than one sensible UMI/detection threshold; include
   cell-type-composition and spatial-geometry-matched nulls.
4. Do not select whichever parameter grid point maximises held-out
   performance. Report failures and uncertainty.

## F. Independent specimen/disease conclusions

For claims contrasting disease states, require >=8 **independent donors**
with >=3 donors in each primary comparison group as an absolute *exploratory*
minimum, and prefer >=10 donors/group or an external independent cohort for
credible biological conclusions. FOVs/sections are NOT independent donors.
The current 4 CosMx-6K specimens CANNOT pass this clinical gate alone.

## Mandatory decision and dataset fallback

- **PASS-EXPLORATORY:** A–C pass; local sheaf and spatial metrics may be
  explored but NOT described as biologically validated.
- **PASS-VALIDATED:** A–E pass, with adequate held-out independent donor
  evidence; F also required before clinical-state interpretation.
- **NO-GO/UNDETERMINED:** Any failed/non-evaluable gate must be reported
  clearly with evidence and constraints. Do not secretly relax thresholds,
  fit on the test set, or call energy minimisation independent prediction.
- Before declaring the CosMx dataset infeasible, check the original paper's
  reference-transferred/protein-gated annotations and linked CosMx Seurat
  object at https://doi.org/10.6084/m9.figshare.31595527, and its source
  processing methods at https://doi.org/10.1038/s41467-026-72482-w.
- If the missing data are structural, search other available datasets with
  true cell-level spatial coordinates, measured gene programmes and matched
  cell identities, >=8 independent donors, sufficiently observed receptor
  pathways, and an independent perturbation/functional outcome. Audit
  access, independent patient count, matched genes, response reference,
  locations and file sizes BEFORE proposing download.
- An especially relevant comparator is 2026 Crohn's Xenium-5K spatial tissue
  (16 tissues from 8 patients) with CellTypist + refined annotations and
  separate stromal cytokine perturbations; GSE307670 / SCP3259,
  https://doi.org/10.1038/s41590-026-02617-0. It is not necessarily
  adequate for a healthy-vs-disease comparison and must be audited first.

## Why these numbers?

- 95% joins is a stringent data-integrity check, not a biological constant.
- 90% broad / 70% subtype annotation prevents an inferred "communication"
  graph dominated by unknown mixtures of cell types.
- >=2 distinct pathways and >=10 signature genes avoids dependence on a
  single sparse or heavily overlapping transcriptional signature.
- >=100 exposed and >=20% of candidate receiver cells exceeds the original
  10/34/8 exposure counts and reduces near-zero-exposure dilution.
- >=5% held-out MSE improvement is a modest but material predeclared
  predictive threshold in noisy data; it is not a universal minimum.
- Downsampling and spatial-radius checks protect against pure
  measurement-density/threshold effects.
- >=8 independent donors is a deliberately *minimal* clinical replication
  gate, not sufficient on its own for publication-strength generalisation.

**Freeze this file before running any new experiments.**
