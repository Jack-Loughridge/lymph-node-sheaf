# Remaining 24 processed recordings

The bounded VM-executor request failed DNS; no repeated downloads were attempted. This does not establish that your normal SSH session cannot reach Dryad. The public processed collection is approximately320MB; the remaining24 total approximately247.15MB by rounded repository sizes. No raw microscopy is required.

First optionally try through normal Windows PowerShell/SSH:

```powershell
ssh tda-vm "cd /home/Jack/Real_Data/lymph_node_sheaf && python3 review/worm_behavior_sheaf_readiness_v1/export_recovery_v3/acquire_remaining_processed.py --download"
```

The helper discovers actual public download links, enforces the exact filename allowlist, caps each object at20MB, checks HDF5 magic/size, and checkpoints private hashes. It fails closed if current public HTML does not expose these links. It has not been verified against a successful network response in this executor. Do not guess file IDs or repeatedly retry a blocked endpoint.

Otherwise open [the public Dryad deposit](https://datadryad.org/dataset/doi:10.5061/dryad.w9ghx3g4v) in your Windows browser. Download the24 names in `REMAINING_FILES.csv`, or download/extract its small processed collection and transfer only these24. Preserve original filenames. The three existing samples and dictionary need no repeat download.

Copy `TRANSFER_REMAINING_WINDOWS.ps1` from this directory to Windows, then run:

```powershell
scp tda-vm:/home/Jack/Real_Data/lymph_node_sheaf/review/worm_behavior_sheaf_readiness_v1/export_recovery_v3/TRANSFER_REMAINING_WINDOWS.ps1 .
powershell -NoProfile -ExecutionPolicy Bypass -File .\TRANSFER_REMAINING_WINDOWS.ps1
```

Enter the Windows download/extraction folder when prompted. The script transfers to the existing private recordings directory, refuses to overwrite existing files, compares Windows/VM SHA-256 after each SCP, and checkpoints a private receipt. If interrupted, preserve that receipt; reconcile already-transferred files before resuming. For an individual allowlisted file the equivalent command is:

```powershell
scp "C:\Downloads\2023-06-29-01-cudata.h5" tda-vm:/home/Jack/Real_Data/lymph_node_sheaf/data/worm_behavior_sheaf_readiness_v1/recordings/
```

After transfer, run `audit_transferred_cohort.py` with the private receipt path using the existing Python environment. This verifies SHA/size, timestamps, actual array bounds and eligible native identities; it writes only private animal/neuron outputs plus aggregate QC. It does not fit a model or certify absent observation masks. Keep all receipts and recordings out of GitHub. Full behavioural and held-out readiness must be reassessed after all actual files pass QC.
