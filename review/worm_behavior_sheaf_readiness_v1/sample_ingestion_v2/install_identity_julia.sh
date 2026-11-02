#!/usr/bin/env bash
# Official, pinned reproduction runtime in ignored project data. No sudo/PATH edits.
set -euo pipefail
worm_root=/home/Jack/Real_Data/lymph_node_sheaf
worm_tools="$worm_root/data/worm_behavior_sheaf_readiness_v1/julia_tools"
worm_archive="$worm_tools/julia-1.9.2-linux-x86_64.tar.gz"
worm_sha=4c2d799f442d7fe718827b19da2bacb72ea041b9ce55f24eee7b1313f57c4383
worm_url=https://julialang-s3.julialang.org/bin/linux/x64/1.9/julia-1.9.2-linux-x86_64.tar.gz
[[ $(uname -m) == x86_64 ]] || { echo 'This pinned binary requires x86_64; stop.' >&2; exit 2; }
[[ -d "$worm_root" ]] || { echo 'Project root absent; stop.' >&2; exit 2; }
mkdir -p "$worm_tools"
[[ ! -L "$worm_tools" && ! -L "$worm_archive" ]] || { echo 'Unexpected symlink; stop.' >&2; exit 2; }
if [[ ! -f "$worm_archive" ]]; then
  echo 'One bounded official-binary request from this normal VM session.'
  [[ ! -e "$worm_archive.part" ]] || { echo 'Existing partial download; inspect instead of overwrite.' >&2; exit 2; }
  if ! curl --fail --location --connect-timeout 8 --max-time 300 --output "$worm_archive.part" "$worm_url"; then
    echo 'Binary download failed. Stop here; use the Windows transfer fallback. No retry.' >&2
    exit 20
  fi
  echo "$worm_sha  $worm_archive.part" | sha256sum --check --status || { echo 'Official checksum mismatch; preserve partial and stop.' >&2; exit 3; }
  mv -- "$worm_archive.part" "$worm_archive"
fi
echo "$worm_sha  $worm_archive" | sha256sum --check
worm_julia="$worm_tools/julia-1.9.2/bin/julia"
if [[ ! -d "$worm_tools/julia-1.9.2" ]]; then
  echo 'Extracting verified portable Julia inside ignored data/.'
  tar --extract --gzip --file "$worm_archive" --directory "$worm_tools" --no-same-owner
fi
[[ -x "$worm_julia" ]] || { echo 'Incomplete Julia directory; inspect instead of overwrite.' >&2; exit 3; }
"$worm_julia" --startup-file=no --version
export JULIA_DEPOT_PATH="$worm_tools/depot"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
mkdir -p "$worm_tools/project"
if "$worm_julia" --startup-file=no --project="$worm_tools/project" -e 'using JLD2; println("JLD2 already usable")'; then
  echo 'Existing isolated JLD2 environment usable.'
else
  echo 'Testing the official Julia package server before package setup.'
  curl --fail --location --connect-timeout 8 --max-time 20 --output "$worm_tools/registries_probe.txt" https://pkg.julialang.org/registries || {
    echo 'Runtime works but package access is blocked. Use the documented Linux-depot fallback; do not decode JLD2 heuristically.' >&2; exit 21;
  }
  timeout 600 "$worm_julia" --startup-file=no --project="$worm_tools/project" -e 'using Pkg; Pkg.add("JLD2"); Pkg.pin("JLD2"); Pkg.status()'
  "$worm_julia" --startup-file=no --project="$worm_tools/project" -e 'using JLD2; println("Native JLD2 ready")'
fi
"$worm_julia" --startup-file=no --project="$worm_tools/project" "$worm_root/review/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2/export_identities_native.jl" "$worm_root/data/worm_behavior_sheaf_readiness_v1/recordings/NeuroPAL_labels.jld2" "$worm_root/data/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2/native_identity"
echo 'Native identity export complete. Raw dictionary and derived labels remain private.'
