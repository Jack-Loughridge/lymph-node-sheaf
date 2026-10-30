#!/usr/bin/env bash
# Public GEO GSE312415 (CosMx) acquisition; no extraction of large archive.
# Modes: probe | download | manifest
set -euo pipefail

GSE="GSE312415"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312415/suppl/GSE312415_files.tar.gz"
ROOT="${IBD_SHEAF_ROOT:-$HOME/Real_Data/spatial_ibd_sheaf}"
ARCHIVES="$ROOT/archives"
METADATA="$ROOT/metadata"
ARCHIVE="$ARCHIVES/${GSE}_files.tar.gz"
MODE="${1:-probe}"
mkdir -p "$ARCHIVES" "$METADATA" "$ROOT/processed" "$ROOT/outputs"

case "$MODE" in
  probe)
    printf 'Source: %s\nDestination: %s\n' "$URL" "$ARCHIVE"
    printf '\nDisk availability:\n'
    df -h "$ROOT"
    printf '\nRemote HTTP response headers (no download):\n'
    curl -fLsSI --connect-timeout 15 --max-time 60 --retry 2 "$URL" \
      | grep -Ei '^(HTTP/|content-length:|content-type:|last-modified:|accept-ranges:)' || true
    printf '\nThis probe does not download the archive.\n'
    ;;
  download)
    command -v wget >/dev/null || { echo 'Error: wget is required.' >&2; exit 1; }
    # GEO reports about 6 GB compressed. Require sufficient room to preserve
    # storage for VM operations; do NOT extract the archive here.
    available=$(df -B1 --output=avail "$ROOT" | tail -n 1 | tr -d '[:space:]')
    existing=0
    if [[ -f "$ARCHIVE" ]]; then existing=$(stat -c %s "$ARCHIVE"); fi
    if (( existing == 0 && available < 18*1024*1024*1024 )); then
      echo "Insufficient free space for a safe initial download (<18 GiB)." >&2
      df -h "$ROOT"
      exit 1
    fi
    if (( available < 10*1024*1024*1024 )); then
      echo 'Insufficient free space to continue safely (<10 GiB free).' >&2
      df -h "$ROOT"
      exit 1
    fi
    printf 'Downloading public processed CosMx GEO archive to %s\n' "$ARCHIVE"
    printf 'The transfer can resume if interrupted. No extraction will occur.\n'
    wget --continue --tries=8 --timeout=90 --waitretry=5 \
      --show-progress --output-document="$ARCHIVE" "$URL"
    printf '\nArchive file size:\n'
    ls -lh "$ARCHIVE"
    printf '\nRemaining disk space:\n'
    df -h "$ROOT"
    ;;
  manifest)
    [[ -s "$ARCHIVE" ]] || { echo 'Archive missing. Run download first.' >&2; exit 1; }
    INDEX="$METADATA/${GSE}_archive_members.txt"
    echo "Reading archive member names (without extracting): $INDEX"
    tar -tzf "$ARCHIVE" > "$INDEX"
    echo "Total entries: $(wc -l < "$INDEX")"
    printf '\nCandidate cell-level expression and spatial-location files:\n'
    (grep -Ei 'exprmat|expression|metadata_file|fov_positions|cell.*(csv|parquet|h5ad)' "$INDEX" || true) | head -n 75
    printf '\nNothing has been extracted. Share the manifest output to choose a small pilot subset.\n'
    ;;
  *)
    echo 'Usage: bash download_ibd_cosmx.sh [probe|download|manifest]' >&2
    exit 2
    ;;
esac
