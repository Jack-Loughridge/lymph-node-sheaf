#!/usr/bin/env python3
"""Summarize Farzad SenNet inventory metadata safely (no downloads)."""
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

def rows(path):
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def size_fmt(n):
    try: n=int(n)
    except (ValueError, TypeError): return "?"
    for unit in ("B","KB","MB","GB","TB"):
        if n < 1024: return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--dir",default="metadata")
    args=ap.parse_args()
    d=Path(args.dir)
    ds=rows(d/"datasets.csv")
    fs=rows(d/"files.csv")
    print(f"Datasets: {len(ds)}; indexed files: {len(fs)}")
    by_type=defaultdict(Counter)
    for r in ds:
        typ=r.get("dataset_type") or "(unknown)"
        access=r.get("access") or "(unknown)"
        by_type[typ][access]+=1
    print("\nDatasets by assay and access:")
    for typ, count in sorted(by_type.items()):
        print(f"  {typ}: "+", ".join(f"{a}: {n}" for a,n in sorted(count.items())))
    files_by_ds=Counter(f.get("sennet_id") for f in fs)
    print("\nIndexed files by dataset (first 12):")
    for sid,n in files_by_ds.most_common(12):
        print(f"  {sid}: {n}")
    print("\nCandidate public data tables (excluding metadata, raw imaging, and channel keys):")
    interesting=[]
    for f in fs:
        path=f.get("rel_path","")
        lower=path.lower()
        if f.get("data_access_level")!="public": continue
        if any(t in lower for t in ("metadata.json", "metadata.tsv", "microscope_",
            ".channels.csv", "/raw/", ".qptiff", "overview", "exposure")):
            continue
        if any(lower.endswith(ext) for ext in (".csv",".tsv",".csv.gz",".tsv.gz",".h5",
           ".h5ad",".mtx",".mtx.gz",".rds",".parquet",".zarr")):
            interesting.append(f)
    for f in sorted(interesting,key=lambda x:x.get("sennet_id","")):
        print("  ",f.get("sennet_id"),size_fmt(f.get("size")),f.get("rel_path"))
    if not interesting:
        print("  None found in public index. The 33 'candidate processed' files may be metadata, not expression matrices.")
    print("\nNote: Protected DBiT records may have additional files omitted from this unauthenticated index.")
    print("This script does not attempt authentication or download data.")

if __name__=="__main__":
    main()
