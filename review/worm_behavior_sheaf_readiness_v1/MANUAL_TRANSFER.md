# Complete the blocked steps using the existing tda-vm alias

## Optional remaining exact archive cleanup

The executor removed the 8,984,620,838-byte Crohn H5AD after checksum/path/open-handle checks. The 6,464,275,219-byte IBD archive could not be removed: the filesystem exposes that directory read-only. It is not an approval rejection or a missing user authorization. No permission escalation or SSH bypass was attempted. A normal VM session may have the required filesystem access.

In Windows PowerShell, after the code is available on the VM:

```powershell
ssh tda-vm "python3 /home/Jack/Real_Data/lymph_node_sheaf/review/worm_behavior_sheaf_readiness_v1/finish_archive_cleanup_on_vm.py"
```

This checks the exact archive, all path components, regular-file/single-link status, size, SHA-256 and visible open users again before unlinking only that file. It leaves code, outputs and archive-member metadata untouched and writes a separate receipt. Do not use a recursive `rm` command. If normal VM process checks reveal a user or status uncertainty, stop. Existing checksum: `32fa2857563a117f974ed507083a996c34a260b52633b30fc3a9e1bc334c7cc6`.

## Browser download: five small processed files

Open the [public Dryad study](https://datadryad.org/dataset/doi:10.5061/dryad.w9ghx3g4v). Save each filename below exactly, into one Windows folder. These links are verified to be the public manifest's destinations; they returned 403 in this web reader. **Successful access in a normal browser has not been verified here.** No fabricated alternative endpoint or repository checksum is supplied.

| Filename | Listed size, rounded | Manifest-linked download |
|---|---:|---|
| `2023-03-30-01-cudata.h5` | 9.71 MB | [Download](https://datadryad.org/downloads/file_stream/4707115) |
| `2023-06-21-01-cudata.h5` | 11.38 MB | [Download](https://datadryad.org/downloads/file_stream/4707119) |
| `2023-06-23-08-cudata.h5` | 10.17 MB | [Download](https://datadryad.org/downloads/file_stream/4707120) |
| `NeuroPAL_labels.jld2` | 41.61 MB | [Download](https://datadryad.org/downloads/file_stream/4707122) |
| `load_data_example.jl` | 1.61 KB | [Download](https://datadryad.org/downloads/file_stream/4707124) |

These three animals were fixed from the manifest before any response inspection. The complete listed collection is about 320 MB, but **do not acquire the rest until this sample passes real-data QC**. The full manifest has 27 HDF5 filenames; the author schema is not a substitute for array inspection.

## Exact transfer commands

Run the supplied script, which uses your configured SSH identity, hostname and port through `tda-vm`:

```powershell
& 'PATH-TO-PROJECT\review\worm_behavior_sheaf_readiness_v1\TRANSFER_WINDOWS.ps1'
```

Or transfer explicitly after checking that the destination does not already contain these filenames:

```powershell
$InputFolder = Read-Host 'Windows folder containing the five files'
$VmDir = '/home/Jack/Real_Data/lymph_node_sheaf/data/worm_behavior_sheaf_readiness_v1/recordings'
ssh tda-vm "mkdir -p '$VmDir'"
Get-FileHash -Algorithm SHA256 "$InputFolder\2023-03-30-01-cudata.h5", "$InputFolder\2023-06-21-01-cudata.h5", "$InputFolder\2023-06-23-08-cudata.h5", "$InputFolder\NeuroPAL_labels.jld2", "$InputFolder\load_data_example.jl"
scp -- "$InputFolder\2023-03-30-01-cudata.h5" "$InputFolder\2023-06-21-01-cudata.h5" "$InputFolder\2023-06-23-08-cudata.h5" "$InputFolder\NeuroPAL_labels.jld2" "$InputFolder\load_data_example.jl" "tda-vm:$VmDir/"
ssh tda-vm "sha256sum '$VmDir/2023-03-30-01-cudata.h5' '$VmDir/2023-06-21-01-cudata.h5' '$VmDir/2023-06-23-08-cudata.h5' '$VmDir/NeuroPAL_labels.jld2' '$VmDir/load_data_example.jl'"
```

Check all five hashes agree. The script automates that comparison and refuses overwrites. A locally computed checksum verifies transfer, not repository authenticity; no trusted repository checksum was available in the accessible listing.

Then run the acquisition script once to validate local formats and record sizes/hashes without redownloading present files:

```powershell
ssh tda-vm "cd /home/Jack/Real_Data/lymph_node_sheaf && python3 review/worm_behavior_sheaf_readiness_v1/acquire_sample.py"
```

The next scientific step is to inspect the actual author loader and JLD2 structure, establish identities and signal/timing provenance, and run the offline numerical audit described in REPRODUCTION.md. File integrity alone does not authorize a model or full acquisition.

If the normal browser also receives 403, contact Dryad support or the authors for these five exact processed exports or refreshed public download access. No messages have been sent on your behalf. This access issue is not a scientific NO-GO.
