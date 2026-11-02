# Small inputs needed to unblock the preregistered analysis

Save original filenames into a Windows folder. These endpoints are the named publication/accession files; successful payload access is not claimed in this executor.

| Input | Exact link |
|---|---|
| S1 authoritative mapping | [41590_2026_2617_MOESM3_ESM.pdf](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41590-026-02617-0/MediaObjects/41590_2026_2617_MOESM3_ESM.pdf) |
| IL1B DEG | [41590_2026_2617_MOESM6_ESM.csv](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41590-026-02617-0/MediaObjects/41590_2026_2617_MOESM6_ESM.csv) |
| IFNG DEG | [41590_2026_2617_MOESM9_ESM.csv](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41590-026-02617-0/MediaObjects/41590_2026_2617_MOESM9_ESM.csv) |
| TNF DEG | [41590_2026_2617_MOESM12_ESM.csv](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41590-026-02617-0/MediaObjects/41590_2026_2617_MOESM12_ESM.csv) |
| Culture replicate matrix | [GSE307526_processed_bulkRNAseq_cytokine_stimulations.csv.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307526/suppl/GSE307526_processed_bulkRNAseq_cytokine_stimulations.csv.gz) |

If a direct link fails, use the named supplements on the [paper](https://www.nature.com/articles/s41590-026-02617-0) or the [GEO supplementary-files listing](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307526). Approximate prior listing sizes: S1 136 KB, each DEG CSV 1.1 MB and culture matrix 1.7 MB; sizes are not independently verified payload measurements.

Run in Windows PowerShell, using the same working address/key used for the atlas:

```powershell
$InputFolder = Read-Host 'Windows folder containing the five supporting files'
$VmTarget = Read-Host 'Working SSH address (Jack@your-VM-address)'
$KeyPath = Read-Host 'Existing SSH private-key path'
$FileNames = @('41590_2026_2617_MOESM3_ESM.pdf','41590_2026_2617_MOESM6_ESM.csv','41590_2026_2617_MOESM9_ESM.csv','41590_2026_2617_MOESM12_ESM.csv','GSE307526_processed_bulkRNAseq_cytokine_stimulations.csv.gz')
ssh -i "$KeyPath" $VmTarget 'mkdir -p /home/Jack/Real_Data/lymph_node_sheaf/data/crohn_xenium/supporting_files'
if ($LASTEXITCODE -ne 0) { throw 'SSH destination setup failed' }
foreach ($Name in $FileNames) {
    $LocalPath = Join-Path $InputFolder $Name
    Get-FileHash -LiteralPath $LocalPath -Algorithm SHA256
    scp -i "$KeyPath" "$LocalPath" "${VmTarget}:/home/Jack/Real_Data/lymph_node_sheaf/data/crohn_xenium/supporting_files/"
    if ($LASTEXITCODE -ne 0) { throw "Transfer failed: $Name" }
}
```

Omit `-i` if using an SSH agent; retain your working non-default port options if applicable. Do not overwrite a different earlier version. Then run the bounded acquisition script: it inspects existing files without redownloading them and records hashes. Compare these with Windows hashes. Presence/integrity does not establish scientific verification.

Additional small exports may require the authors: run-specific gene_panel.json or analysis_summary.html; negative-control probe rows with specimen/cell indexing and target names; independent section outlines and anatomical labels registered to centroids. Do not transfer whole raw imaging or the full transcript cloud merely to recover these controls. A concise request would ask for these exact derived exports, clarify whether receptor competence is demonstrated in primary intestinal fibroblasts for the three individual cytokines, and request the culture replicate treatment/control metadata. No author message has been sent.
