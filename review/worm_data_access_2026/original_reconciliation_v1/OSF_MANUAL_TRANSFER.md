# Original processed exports: concise manual transfer

One fresh attempt per OSF link failed DNS in the execution sandbox; the web reader also could not retrieve them. This does not establish that Windows or a normal VM SSH session is blocked. No more endpoint retries or nine-animal downloads were attempted.

In Windows PowerShell, inspect the final response headers/file size first (or inspect the OSF file page if Content-Length is absent). These are the processed exports, not raw microscopy:

```powershell
curl.exe --head --location "https://osf.io/9mecf/download"
curl.exe --head --location "https://osf.io/ty6sp/download"
```

After checking sizes, download to a new directory and verify the **reported version-1** checksums:

```powershell
$ExportFolder = Join-Path $env:USERPROFILE "Downloads\WormOriginalExports"
New-Item -ItemType Directory -Force -Path $ExportFolder | Out-Null
curl.exe --fail --location "https://osf.io/9mecf/download" --output "$ExportFolder\exported_data.tar.gz"
if ($LASTEXITCODE -ne 0) { throw "WT download failed" }
curl.exe --fail --location "https://osf.io/ty6sp/download" --output "$ExportFolder\exported_data_unc31.tar.gz"
if ($LASTEXITCODE -ne 0) { throw "unc-31 download failed" }
if ((Get-FileHash "$ExportFolder\exported_data.tar.gz" -Algorithm MD5).Hash.ToLowerInvariant() -ne "7ae518fea414f664e2493238afe209bb") { throw "WT version/checksum differs; investigate before use" }
if ((Get-FileHash "$ExportFolder\exported_data_unc31.tar.gz" -Algorithm MD5).Hash.ToLowerInvariant() -ne "e5c83242baa131bca95ee58bc7b58bec") { throw "unc-31 version/checksum differs; investigate before use" }

$VmTarget = "Jack@YOUR_VM_PUBLIC_IP_OR_SSH_DNS"
$KeyPath = "C:\Users\YOUR_WINDOWS_USERNAME\.ssh\YOUR_AZURE_KEY.pem"
ssh -i "$KeyPath" $VmTarget "mkdir -p /home/Jack/Real_Data/lymph_node_sheaf/data/worm_original_exports"
scp -i "$KeyPath" "$ExportFolder\exported_data.tar.gz" "${VmTarget}:/home/Jack/Real_Data/lymph_node_sheaf/data/worm_original_exports/"
scp -i "$KeyPath" "$ExportFolder\exported_data_unc31.tar.gz" "${VmTarget}:/home/Jack/Real_Data/lymph_node_sheaf/data/worm_original_exports/"
```

Use the endpoint/key that successfully transferred the NWBs. Omit `-i` if using your existing SSH agent. For a nondefault SSH port use `-p PORT` for ssh and `-P PORT` for scp. Windows key and public VM address are not known in this session. Use a new directory if archives already exist rather than overwriting an earlier version. MD5 values are from Encore's reproduction record, not hashes measured from original downloads here; the inspection script additionally records SHA-256 locally. No credentials or tokens are required for the advertised public OSF links.

On the VM, reconcile the WT flat export without extracting files or loading pickle objects:

```bash
cd /home/Jack/Real_Data/lymph_node_sheaf
/home/Jack/micromamba/envs/humann_py312/bin/python \
  review/worm_data_access_2026/original_reconciliation_v1/reconcile.py \
  --wt-archive data/worm_original_exports/exported_data.tar.gz
```

The scanner matches `_ds_name.txt` acquisition names to **20220331_114306, 20230309_143050 and 20230314_164209**. The numbered OSF export prefix is not a DANDI animal ID. It validates dimensions, original stimulus-volume alignment and unique labels, and retains mismatches rather than silently shifting indices. It writes event-level observations only to ignored `data/worm_data_access_2026/original_reconciliation_v1/`.

The flat text exports do not, by their documented writer, include original NaN masks, matchless provenance, response-list flags or all manual corrections. Dense `_gcamp.txt` values must not automatically be labelled measured observations. The complete `exported_data_full.tar.gz` is a distinct library-compatible export; its direct GUID and size remain unverified. An original `fconn.pickle` may contain inclusion state, but our scanner deliberately does not execute it. A safe author-exported JSON of the relevant attributes is preferable.
