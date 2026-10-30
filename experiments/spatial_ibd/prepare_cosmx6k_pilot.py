#!/usr/bin/env python3
"""Stream only the 6K CosMx pilot files from GEO GSE312415.

No full-archive extraction. Expression matrices are recompressed as .csv.gz
on the fly so the VM never holds the large uncompressed matrices.
Designed for public data; retain all data on the VM, outside Git.
"""
import argparse
import csv
import gzip
import os
from pathlib import Path
import shutil
import tarfile

PREFIXES = (
    "Cos1400R355RRNA6K4u1716",  # healthy + UC noninflamed
    "Cos1380R387RRNA6K4u2716",  # two UC inflamed specimens
)
SUFFIXES = ("_exprMat_file.csv", "_metadata_file.csv", "_fov_positions_file.csv")
SAMPLES = [
    ("GSM9345416", "Rectum_HC", "Cos1400R355RRNA6K4u1716", 32, 65),
    ("GSM9345417", "Rectum_UC_Noninflamed", "Cos1400R355RRNA6K4u1716", 1, 31),
    ("GSM9345418", "Rectum_UC_Inflamed", "Cos1380R387RRNA6K4u2716", 18, 45),
    ("GSM9345419", "Rectum_UC_Inflamed", "Cos1380R387RRNA6K4u2716", 1, 17),
]
MIB = 1024 * 1024
GIB = 1024 ** 3

def inspect(root):
    for prefix in PREFIXES:
        print(f"\n{prefix}:")
        for suffix in SUFFIXES:
            p = root / "counts" / (prefix + suffix + ".gz") if suffix.startswith("_expr") else root / "metadata" / (prefix + suffix)
            if not p.is_file():
                print(f"  Missing: {p.name}")
                continue
            with gzip.open(p, "rt", newline="", encoding="utf-8-sig") if suffix.startswith("_expr") else p.open("rt", newline="", encoding="utf-8-sig") as fh:
                reader = csv.reader(fh)
                header = next(reader, [])
            print(f"  {p.name}: {p.stat().st_size / MIB:,.1f} MiB; "
                  f"{len(header)} columns; first columns: {header[:12]}")
    print("\nGEO specimen-to-FOV mapping is in metadata/cosmx6k_samples.csv.")
    print("Be careful: the two expression files contain FOUR different specimens.")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive", type=Path, default=Path.home() / "Real_Data/spatial_ibd_sheaf/archives/GSE312415_files.tar.gz")
    p.add_argument("--out", type=Path, default=Path.home() / "Real_Data/spatial_ibd_sheaf/pilot_6k")
    p.add_argument("--with-expression", action="store_true", help="Also stream/recompress the two 6K expression matrices (more disk and CPU)")
    p.add_argument("--inspect", action="store_true", help="Inspect already extracted CSV headers without reading archive")
    p.add_argument("--reserve-gib", type=float, default=10, help="Minimum free disk to preserve (default 10 GiB)")
    args=p.parse_args()
    root=args.out.expanduser()
    (root/"metadata").mkdir(parents=True,exist_ok=True)
    (root/"counts").mkdir(parents=True,exist_ok=True)
    metadata=root/"metadata/cosmx6k_samples.csv"
    with metadata.open("w",newline="") as fh:
        writer=csv.writer(fh)
        writer.writerow(("geo_sample","condition","archive_prefix","fov_start","fov_end"))
        writer.writerows(SAMPLES)
    if args.inspect:
        inspect(root)
        return
    archive=args.archive.expanduser()
    if not archive.is_file():
        raise SystemExit(f"Cannot find archive: {archive}")
    free=shutil.disk_usage(root).free
    if free < args.reserve_gib*GIB+2*GIB:
        raise SystemExit(f"Only {free/GIB:.1f} GiB free; need >= {args.reserve_gib+2:.1f} GiB before extraction.")
    want={prefix+suffix for prefix in PREFIXES for suffix in SUFFIXES
          if args.with_expression or not suffix.startswith("_expr")}
    found=set()
    print(f"Streaming from: {archive}",flush=True)
    print("Will extract:",*sorted(want),sep="\n  ",flush=True)
    # Sequential gzip tar streaming: does not unpack unrelated samples.
    with tarfile.open(archive,"r|gz") as tf:
        for member in tf:
            name=Path(member.name).name
            if name not in want or not member.isfile():
                continue
            if name in found:
                raise RuntimeError(f"Duplicate member {name}; refusing to overwrite.")
            found.add(name)
            isexpr=name.endswith("_exprMat_file.csv")
            dest=(root/"counts"/(name+".gz")) if isexpr else (root/"metadata"/name)
            if dest.exists():
                print(f"Already present, skipping: {dest.name}",flush=True)
                continue
            if member.size > 30*GIB:
                raise RuntimeError(f"Unusually large file {name}: {member.size/GIB:.1f} GiB")
            print(f"Extracting {name}; original {member.size/GIB:.2f} GiB, "
                  f"{'gzip-compressed' if isexpr else 'unchanged'}...",flush=True)
            temp=dest.with_name(dest.name+".part")
            input_file=tf.extractfile(member)
            if input_file is None:
                raise RuntimeError("Archive member could not be read: "+name)
            try:
                with temp.open("wb") as raw_out:
                    if isexpr:
                        dst=gzip.GzipFile(fileobj=raw_out,mode="wb",compresslevel=5,mtime=0)
                    else:
                        dst=raw_out
                    try:
                        while True:
                            chunk=input_file.read(4*MIB)
                            if not chunk: break
                            dst.write(chunk)
                            if shutil.disk_usage(root).free < args.reserve_gib*GIB:
                                raise RuntimeError(f"Stopping: free space fell below {args.reserve_gib} GiB")
                    finally:
                        if isexpr:
                            dst.close()
                temp.rename(dest)
                print(f"  Saved {dest.name}: {dest.stat().st_size/GIB:.2f} GiB",flush=True)
            except BaseException:
                temp.unlink(missing_ok=True)
                raise
    missing=want-found
    if missing:
        raise RuntimeError("Expected files absent from archive: "+", ".join(sorted(missing)))
    print("\nSuccess: selected files extracted. Disk available:",
          f"{shutil.disk_usage(root).free/GIB:.1f} GiB",flush=True)
    inspect(root)

if __name__=="__main__":
    main()
