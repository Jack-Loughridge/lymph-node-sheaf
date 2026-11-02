# One minimal manual transfer to resolve the strongest public-file uncertainty

No full12GB/68GB data download is recommended. Executor requests to Dryad, Springer media and Yale failed DNS once per service; the web reader accesses manifests but cannot deliver the Ruhl binary ZIP. This does not prove a failure of the normal Azure VM or absence of data. The GitHub connector can read author code, which is why a code-level export-schema audit was possible; it is not a substitute for physiological files.

Public Source Data link from the original paper:

https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-76758-z/MediaObjects/41467_2026_76758_MOESM10_ESM.zip

The actual file size and repository checksum are NOT VERIFIED. In Windows, use the paper's **Source Data** link or the exact URL above. Download only if the browser indicates a reasonably small file; the initial audit accepts≤50,000,000bytes. If larger, report the displayed size/manifest first; do not transfer a huge archive. Do not download raw movies or all supplementary material. Preserve the original filename.

After downloading, Windows PowerShell:

```powershell
$File = Join-Path $HOME 'Downloads\41467_2026_76758_MOESM10_ESM.zip'
if ((Get-Item -LiteralPath $File).Length -gt 50000000) { throw 'Stop: exceeds initial 50MB budget' }
$Hash = (Get-FileHash -LiteralPath $File -Algorithm SHA256).Hash.ToLowerInvariant()
$Receipt = @{algorithm='SHA256';filename=(Split-Path $File -Leaf);bytes=(Get-Item -LiteralPath $File).Length;sha256=$Hash}
$ReceiptPath = Join-Path $HOME 'Downloads\ruhl_transfer_PRIVATE.json'
if (Test-Path -LiteralPath $ReceiptPath) { throw 'Preserve existing receipt; choose a versioned name' }
$Receipt | ConvertTo-Json | Set-Content -Encoding utf8 $ReceiptPath
ssh tda-vm 'test ! -e /home/Jack/Real_Data/lymph_node_sheaf/data/perturbational_bioelectric_integration_2026/41467_2026_76758_MOESM10_ESM.zip'
if ($LASTEXITCODE -ne 0) { throw 'Destination already exists; inspect rather than overwrite' }
scp "$File" tda-vm:/home/Jack/Real_Data/lymph_node_sheaf/data/perturbational_bioelectric_integration_2026/
if ($LASTEXITCODE -ne 0) { throw 'SCP failed' }
$Remote = ssh tda-vm 'sha256sum /home/Jack/Real_Data/lymph_node_sheaf/data/perturbational_bioelectric_integration_2026/41467_2026_76758_MOESM10_ESM.zip'
if ($LASTEXITCODE -ne 0 -or ($Remote -split '\s+')[0] -ne $Hash) { throw 'Transfer hash mismatch' }
scp "$ReceiptPath" tda-vm:/home/Jack/Real_Data/lymph_node_sheaf/data/perturbational_bioelectric_integration_2026/
```

The receipt is a sender/receiver transfer check, not a repository-provided hash. After it arrives, verify bytes/hash from that private receipt and run:

```bash
cd /home/Jack/Real_Data/lymph_node_sheaf
PY=/home/Jack/micromamba/envs/humann_py312/bin/python
"$PY" review/perturbational_bioelectric_integration_2026/audit_small_payload.py data/perturbational_bioelectric_integration_2026/41467_2026_76758_MOESM10_ESM.zip --sha256 WINDOWS_HASH --tag ruhl_v1
```

The auditor inspects ZIP CRC/paths and bounded XLSX/XML, MAT/HDF5 or text schemas without extracting everything or executing archive code. Igor/binary formats remain explicitly undecoded if no verified reader is installed. It never emits automatic DATA GO. Trace identity, simultaneous grouping, culture units, spatial coordinates, conditions and input logs still require source-specific reconciliation. If the file is only figure aggregates or selected ungrouped traces, stop the full regional-model route and request the missing processed fields. Do not turn an exemplar into an independent cohort.

No need to download Cahill's62MB events ZIP merely to rediscover the continuous-state/full-roster omission already demonstrated in its exporter. The alternative minimal action is a representative author trace/ROI export described in `AUTHOR_REQUESTS.md`.
