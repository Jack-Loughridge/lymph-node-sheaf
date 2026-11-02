# Remaining small external inputs

The complete project data/output search found no local GSE307526 count CSV, cytokine DEG table or Supplementary Table S1. Fresh requests failed DNS in the execution sandbox; publisher/NCBI web retrieval also failed. No full spatial data were redownloaded and no repeated download loop was started. These failures do not establish that a normal Windows browser or VM SSH terminal cannot retrieve the files.

Use the following actual accession/supplement files. Approximate sizes are from the publisher/GEO listings, not payload sizes independently measured here.

| Input | File | Approximate size | Download |
|---|---|---:|---|
| Culture replicate expression | `GSE307526_processed_bulkRNAseq_cytokine_stimulations.csv.gz` | 1.7 MB | [NCBI](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307526/suppl/GSE307526_processed_bulkRNAseq_cytokine_stimulations.csv.gz) |
| Authoritative patient/tissue cohort table, Supplementary Table 1 | `41590_2026_2617_MOESM3_ESM.pdf` | 135.5 KB | [Publisher](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41590-026-02617-0/MediaObjects/41590_2026_2617_MOESM3_ESM.pdf) |
| IL1B DESeq2 results, Supplementary Table 4 | `41590_2026_2617_MOESM6_ESM.csv` | 1.1 MB | [Publisher](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41590-026-02617-0/MediaObjects/41590_2026_2617_MOESM6_ESM.csv) |
| IFNG DESeq2 results, Supplementary Table 7 | `41590_2026_2617_MOESM9_ESM.csv` | 1.1 MB | [Publisher](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41590-026-02617-0/MediaObjects/41590_2026_2617_MOESM9_ESM.csv) |
| TNF DESeq2 results, Supplementary Table 10 | `41590_2026_2617_MOESM12_ESM.csv` | 1.1 MB | [Publisher](https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41590-026-02617-0/MediaObjects/41590_2026_2617_MOESM12_ESM.csv) |

Primary listings: [GSE307526](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307526), [paper and supplements](https://www.nature.com/articles/s41590-026-02617-0). If direct supplement links fail in a browser, use the named Supplementary Table links on the paper page. The corresponding PMC media filenames were resolved from the primary supplement links; no successful payload download is claimed.

Save these files in a dedicated local folder, leaving the GEO file gzipped. Preserve source filenames. Copy only these small supporting files using the same SSH connection already used for the atlas, for example in Windows PowerShell:

```powershell
$Support = Read-Host "Full Windows folder containing these supporting files"
$VmTarget = Read-Host "Working SSH address, e.g. Jack@your-VM-address"
$KeyPath = Read-Host "Your existing private-key filepath"
ssh -i "$KeyPath" $VmTarget "mkdir -p /home/Jack/Real_Data/lymph_node_sheaf/data/crohn_xenium/supporting_files"
foreach ($File in Get-ChildItem -LiteralPath $Support -File) {
    Get-FileHash -LiteralPath $File.FullName -Algorithm SHA256
    scp -i "$KeyPath" "$($File.FullName)" "${VmTarget}:/home/Jack/Real_Data/lymph_node_sheaf/data/crohn_xenium/supporting_files/"
    if ($LASTEXITCODE -ne 0) { throw "Supporting-file transfer failed" }
}
```

Omit `-i` for an SSH agent, or use the working port settings. Do not overwrite an earlier version. Keep S1, workbook-derived mappings and all individual culture/spatial data under ignored `data/`; none is part of the public report.

The count CSV is needed to examine scale, actual treatment/control column identities and replicate reproducibility. The DESeq2 result tables may provide signed effects and adjusted p-values but cannot independently establish replicate reproducibility. Four culture replicates are not four spatial patients, and a stromal cell-line dictionary is not an organoid donor-response validation.

These files will complete provenance and response-map assessment. They cannot, by themselves, repair sparse receptor-subunit co-detection already measured in the spatial object. Predictive fitting remains blocked by both unresolved inputs and the measured coverage failure.
