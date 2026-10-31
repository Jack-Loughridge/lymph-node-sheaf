#!/usr/bin/env python3
"""Audit CosMx 6K pilot cell counts, coordinates and illustrative ligand/receptor coverage.

Reads the two metadata CSVs and only the headers of compressed expression
matrices. Requires only Python's standard library. No download, no modelling.
"""
import argparse
import collections
import csv
import gzip
import json
from pathlib import Path

PREFIXES = {
    "Cos1400R355RRNA6K4u1716": [
        ("GSM9345416", "Rectum_HC", 32, 65),
        ("GSM9345417", "Rectum_UC_Noninflamed", 1, 31),
    ],
    "Cos1380R387RRNA6K4u2716": [
        ("GSM9345418", "Rectum_UC_Inflamed", 18, 45),
        ("GSM9345419", "Rectum_UC_Inflamed", 1, 17),
    ],
}
# Examples ONLY. Replace with a curated, version-pinned human LR database
# before fitting any communication model. Detection/expression not yet checked.
PAIRS = [
    ("CXCL12","CXCR4","Chemokine"),
    ("CXCL13","CXCR5","Chemokine"),
    ("CCL19","CCR7","Chemokine"),
    ("CCL20","CCR6","Chemokine"),
    ("CXCL8","CXCR1","Chemokine"),
    ("CXCL8","CXCR2","Chemokine"),
    ("IL6","IL6R","Cytokine"),
    ("IL1B","IL1R1","Cytokine"),
    ("TNF","TNFRSF1A","Cytokine"),
    ("IFNG","IFNGR1","Cytokine"),
    ("TGFB1","TGFBR2","Growth factor"),
    ("OSM","OSMR","Cytokine"),
    ("EGF","EGFR","Growth factor"),
    ("JAG1","NOTCH1","Contact"),
    ("DLL4","NOTCH1","Contact"),
]

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=Path,
        default=Path.home()/"Real_Data"/"spatial_ibd_sheaf"/"pilot_6k")
    args=ap.parse_args()
    root=args.root.expanduser()
    report={"samples":[],"gene_panels":{},"coordinate_headers":{}}
    for prefix,specimens in PREFIXES.items():
        matrix=root/"counts"/(prefix+"_exprMat_file.csv.gz")
        metadata=root/"metadata"/(prefix+"_metadata_file.csv")
        if not matrix.exists() or not metadata.exists():
            raise SystemExit(f"Missing {matrix} or {metadata}")
        with gzip.open(matrix,"rt",encoding="utf-8-sig",newline="") as fh:
            hdr=next(csv.reader(fh))
        genes={g.upper() for g in hdr[3:]}
        report["gene_panels"][prefix]={
            "columns":len(hdr), "gene_columns":len(genes),
            "example_pairs_with_both_genes":[list(p) for p in PAIRS
                                           if p[0] in genes and p[1] in genes],
            "example_pairs_missing_genes":[list(p) for p in PAIRS
                                           if p[0] not in genes or p[1] not in genes],
        }
        print(f"\n{prefix}: {len(genes)} gene columns")
        print(f"  Example LR pairs covered (presence in panel only): "
              f"{len(report['gene_panels'][prefix]['example_pairs_with_both_genes'])}/{len(PAIRS)}")
        with metadata.open(newline="",encoding="utf-8-sig") as fh:
            reader=csv.DictReader(fh)
            fields=reader.fieldnames or []
            print("  Metadata columns:", ", ".join(fields))
            coordinate_fields=[k for k in fields if
                "center" in k.lower() or "global" in k.lower()
                or k.lower() in ("fov","cell_id","cell","slide")]
            report["coordinate_headers"][prefix]=coordinate_fields
            fov_key=next((k for k in fields if k.lower()=="fov"),None)
            if not fov_key:
                raise RuntimeError(f"No FOV column in {metadata}")
            fov_counts=collections.Counter()
            missing_fov=0
            for row in reader:
                try: fov=int(float(row[fov_key]))
                except (ValueError,TypeError):
                    missing_fov+=1
                    continue
                fov_counts[fov]+=1
        print(f"  Metadata cells: {sum(fov_counts.values()):,}, FOVs: {len(fov_counts)}; "
              f"unparsed FOV: {missing_fov}")
        print("  Coordinate columns:",coordinate_fields)
        if fov_counts:
            print(f"  FOV range: {min(fov_counts)}–{max(fov_counts)}")
        for geo,condition,lo,hi in specimens:
            count=sum(n for f,n in fov_counts.items() if lo<=f<=hi)
            rec={"geo":geo,"condition":condition,"prefix":prefix,
                 "fov_start":lo,"fov_end":hi,"metadata_cell_count":count}
            report["samples"].append(rec)
            print(f"    {geo} {condition}, FOV {lo}–{hi}: {count:,} cells")
        uncovered={str(f):n for f,n in sorted(fov_counts.items())
                   if not any(lo<=f<=hi for _,_,lo,hi in specimens)}
        if uncovered:
            print("  WARNING: FOVs without a GEO specimen mapping:",uncovered)
            report.setdefault("unmapped_fovs",{})[prefix]=uncovered
    output=root/"metadata"/"initial_audit.json"
    output.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(f"\nSaved {output}")
    print("IMPORTANT: panel coverage is not expression detection. "
          "Do not infer signalling or spatial distances from this audit.")

if __name__=="__main__":
    main()
