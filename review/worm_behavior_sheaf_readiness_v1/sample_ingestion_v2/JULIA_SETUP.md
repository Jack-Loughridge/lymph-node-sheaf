> Recovery update: Julia 1.9.2/JLD2 0.6.7 are installed and native validation now succeeds. Do not reinstall. See [recovery reproduction](../export_recovery_v3/REPRODUCTION.md) and the consolidated report. The setup instructions below are historical.

# Native Julia identity decoding: minimal isolated setup

All five Windows/VM hashes and sizes match. Do not repeat the dataset downloads.

## What could be tested here

No Julia executable, juliacall runtime or cached Julia distribution was found. One official-runtime metadata request failed DNS. This executor does **not** have your `tda-vm` SSH alias configured, so it cannot execute your normal Windows SSH session or establish whether that session has public-network access. This does not establish VM-wide network failure. No credential configuration, proxy bypass or Python environment change was attempted.

The VM architecture visible here is **x86_64**. Use the official **Linux glibc x86_64** archive, not a Windows installer. Julia **1.9.2** matches the author's documented loading environment; it is an older, unmaintained reproduction runtime, isolated for this audit rather than a general system installation. The [official release listing](https://julialang.org/downloads/oldreleases/) supplies the checksum below. [Official binary guidance](https://julialang.org/downloads/platform/) permits a portable directory install without sudo or PATH edits.

## Preferred: one command in normal Windows PowerShell

```powershell
ssh tda-vm "bash /home/Jack/Real_Data/lymph_node_sheaf/review/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2/install_identity_julia.sh"
```

This tries the official binary download once, verifies SHA-256, extracts under ignored project data, tests the package server, installs and pins **only JLD2** in an isolated project/depot, then runs native identity validation. SHA/TOML are Julia standard libraries. HDF5.jl, GPU/imaging libraries, PyCall and reconstruction packages are unnecessary for the identity-loading block. Existing Python environments and shell configuration remain unchanged.

All runtime/package files stay under:

`/home/Jack/Real_Data/lymph_node_sheaf/data/worm_behavior_sheaf_readiness_v1/julia_tools/`

A failed binary request stops with exit 20; a failed package-server probe stops with exit 21. No repeated endpoint loop or automatic Julia-version upgrade occurs. If the native metadata validation fails, retain the actual error for reconciliation; do not replace it with Python object decoding.

## Windows download and SCP fallback for the Linux runtime

If the normal SSH download fails, run in PowerShell:

```powershell
$DownloadFolder = Read-Host 'Windows folder for the Julia Linux archive'
$JuliaArchive = Join-Path $DownloadFolder 'julia-1.9.2-linux-x86_64.tar.gz'
$ExpectedHash = '4c2d799f442d7fe718827b19da2bacb72ea041b9ce55f24eee7b1313f57c4383'
if (Test-Path -LiteralPath $JuliaArchive) { throw 'Archive exists; verify it instead of overwriting' }
Invoke-WebRequest -Uri 'https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.2-linux-x86_64.tar.gz' -OutFile $JuliaArchive
$ActualHash = (Get-FileHash -LiteralPath $JuliaArchive -Algorithm SHA256).Hash.ToLowerInvariant()
if ($ActualHash -ne $ExpectedHash) { throw 'Official archive checksum mismatch' }
$VmTools = '/home/Jack/Real_Data/lymph_node_sheaf/data/worm_behavior_sheaf_readiness_v1/julia_tools'
ssh tda-vm "mkdir -p '$VmTools'; test ! -e '$VmTools/julia-1.9.2-linux-x86_64.tar.gz'"
if ($LASTEXITCODE -ne 0) { throw 'Destination exists or cannot be prepared; inspect first' }
scp -- "$JuliaArchive" "tda-vm:$VmTools/"
if ($LASTEXITCODE -ne 0) { throw 'Archive transfer failed' }
ssh tda-vm "sha256sum '$VmTools/julia-1.9.2-linux-x86_64.tar.gz'"
ssh tda-vm "bash /home/Jack/Real_Data/lymph_node_sheaf/review/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2/install_identity_julia.sh"
```

The script rechecks the official digest on the VM, then installs from the transferred archive. If an archive already exists, hash-check it rather than repeating the download. No dataset copies are overwritten.

## If the runtime works but Julia package downloads also fail

Copying the runtime alone does not provide JLD2. If **an existing WSL Linux installation** is available on your Windows machine, prepare the same Linux runtime and project/depot there, with working network access. Do not install WSL solely for this small audit without assessing the extra download/storage cost.

Within WSL, with the verified Linux archive in a staging directory:

```bash
mkdir -p "$HOME/worm_jld2_stage"
tar -xzf /mnt/c/YOUR_DOWNLOAD_FOLDER/julia-1.9.2-linux-x86_64.tar.gz -C "$HOME/worm_jld2_stage"
export JULIA_DEPOT_PATH="$HOME/worm_jld2_stage/depot"
export JULIA_PKG_PRECOMPILE_AUTO=0
mkdir -p "$HOME/worm_jld2_stage/project"
"$HOME/worm_jld2_stage/julia-1.9.2/bin/julia" --startup-file=no --project="$HOME/worm_jld2_stage/project" -e 'using Pkg; Pkg.add("JLD2"); Pkg.pin("JLD2"); using JLD2; Pkg.status()'
tar -czf /mnt/c/YOUR_DOWNLOAD_FOLDER/worm_jld2_linux_depot.tar.gz -C "$HOME/worm_jld2_stage" project depot
```

This is a **Linux** package depot, including any required Linux artifacts; do not copy Windows compiled packages as a Linux installation. Compute its SHA-256 in PowerShell, SCP it via `tda-vm`, compare the VM hash, and inspect the two archive roots before extracting into a fresh/empty `julia_tools` package directory. Preserve an earlier failed depot rather than overwriting it. Re-run the installer; its native `using JLD2` check avoids online installation once the depot is usable. If no Linux staging environment is available, report the package-server error; native decoding remains blocked rather than switching to a heuristic reader.

## Native validation and private output

`export_identities_native.jl` uses the actual transferred author loading operation:

```julia
f = jldopen(dict_path, "r")
NeuroPAL_labels = f["dict_match_dict"]
close(f)
```

It checks the trace-row lookup against the neuron-class lookup for label, class, confidence and laterality; checks integer indices and conflicts; then writes a private TSV and a TOML receipt containing Julia/JLD2 versions and hashes. It does not convert class names into individual identities or use the HDF5 registration score as the manual NeuroPAL confidence.

Outputs are under ignored data:

- `sample_ingestion_v2/native_identity/labels_PRIVATE.tsv`
- `sample_ingestion_v2/native_identity/native_receipt_PRIVATE.toml`

After these exist, resume the prepared numerical/identity audit. No full sheaf or higher-level compatibility metric should be fitted by this setup procedure. The full 27-animal acquisition remains conditional on completion of the three-sample technical/provenance checks.
