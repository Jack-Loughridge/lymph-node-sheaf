# Background method audit — original results preserved

**Methodological discrepancy confirmed.** The frozen screen is a conservative research rule, not the 10x per-control background estimator. It pools each cell's negative-control probe and codeword totals, then requires a gene count >=max(2,1+tissue pooled-control 99.9th percentile) and greater than that cell's pooled controls. It therefore combines different control mechanisms, omits the number of control targets and imposes a floor unrelated to per-gene background expectation. Its high percentile and minimum-two requirement can censor real low counts. It was labelled an aggregate-control screen in the original report, but must not be described as 10x-recommended gene positivity.

The official metric instead divides negative-control **probe** counts by cells and control targets. Multiplying this mean by target-gene count estimates an average false-positive burden; it does not supply a universal cell/gene expression cutoff. Codeword controls monitor decoding and are reported separately. The guidance describes a run-level QC metric, not a requirement to discard all one-molecule ligand detections. [10x analysis-summary documentation](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/1.9/analysis/xoa-output-analysis-summary), [control feature definitions](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html).

Let P_s be the sum of cell-assigned probe-control counts, N_s cells, K_s verified negative-control probe targets and G_s target genes. Then b_s=P_s/(N_s K_s) and estimated false positives per cell=b_s G_s. Summed counts plus a verified K_s suffice for the global mean, **but not** per-control heterogeneity, regional background or gene-specific false-positive rates. Probe-target count is not the number of blank codewords or biological genes.

## Actual data audit

The exported H5AD retains separate per-cell control_probe_counts and control_codeword_counts, but zero individually named negative-control features. Across all 2,131,098 cells these sum to 6,745 probe and 5,197 codeword counts. Across the frozen 1,425,592 QC cells: 5,379 probe and 4,466 codeword counts. No source data or old exposure calls were modified.

As an explicitly **conditional** calculation using 20 control targets from general 10x documentation, tissue means range 0.00005210–0.00030947 probe counts/control/cell (all exported cells), or 0.00006487–0.00035695 in QC cells. Multiplying by the exported 5,001 genes yields 0.2605–1.5477 and 0.3244–1.7851 estimated false-positive transcripts/cell respectively. This is not a verified run-specific estimate: the true panel/codebook denominator and original target universe are missing, and the filtered export is not necessarily the whole original run. Low aggregate controls do not exclude individual off-target probes, transcript misassignment or dropout. All values are reproduced in background_metrics.csv with conditional flags.

Obtain small run-specific gene_panel.json/analysis_summary.html and negative-control feature counts before making calibrated background claims. Do not assume the typical 20-control denominator is verified. If a statistical per-neighbourhood test is used, declare it as this study's model and validate against control spatial distributions; do not attribute it to 10x.

## Effect on interpretation of the frozen NO-GO

No frozen result is recomputed or reclassified. Its complete-receptor >=1-molecule sensitivity removed the disputed threshold altogether; IFNG/IL1B coverage still failed in audited populations. Thus the screen discrepancy is real, but correcting background alone cannot repair the original measured-complex coverage result. This does not exclude protein-mediated signalling without detectable receptor transcripts.

The proposed context representation avoids per-cell receptor detection on biological grounds and is a distinct question, not a relaxed rerun of the original gate. Prediction of an external response must be demonstrated independently before any compatibility model is accepted.
