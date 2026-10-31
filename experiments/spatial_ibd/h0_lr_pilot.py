#!/usr/bin/env python3
"""Exploratory 0D persistence of spatial ligand/receptor *potential* in CosMx 6K.

This is NOT a sheaf compatibility score or COMMOT result. Edges are present
when sender ligand & receiver receptor are both detected within a spatial radius.
Filtration strength is the number of supported LR pairs, never raw expression.
One representative FOV per GEO specimen. Does not mix separate specimens.
"""
import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

import numpy as np
from scipy.spatial import cKDTree

PREFIXES = {
    "Cos1400R355RRNA6K4u1716": [
        ("GSM9345416", "healthy", 32, 65),
        ("GSM9345417", "uc_noninflamed", 1, 31),
    ],
    "Cos1380R387RRNA6K4u2716": [
        ("GSM9345418", "uc_inflamed", 18, 45),
        ("GSM9345419", "uc_inflamed", 1, 17),
    ],
}

# DIAGNOSTIC ONLY: replace with curated, versioned LR pairs before biology.
# IL6/IFNG/TGFB signalling, for example, generally involve additional
# receptor subunits not represented by these simplistic pair entries.
PAIRS = [
    ("CXCL12", "CXCR4", "chemokine"),
    ("CXCL13", "CXCR5", "chemokine"),
    ("CCL19", "CCR7", "chemokine"),
    ("CCL20", "CCR6", "chemokine"),
    ("CXCL8", "CXCR1", "chemokine"),
    ("CXCL8", "CXCR2", "chemokine"),
    ("IL6", "IL6R", "cytokine"),
    ("IL1B", "IL1R1", "cytokine"),
    ("TNF", "TNFRSF1A", "cytokine"),
    ("IFNG", "IFNGR1", "cytokine"),
    ("TGFB1", "TGFBR2", "growth_factor"),
    ("OSM", "OSMR", "cytokine"),
    ("EGF", "EGFR", "growth_factor"),
    ("JAG1", "NOTCH1", "contact"),
    ("DLL4", "NOTCH1", "contact"),
]


class DSU:
    def __init__(self, n):
        self.parent = np.arange(n, dtype=np.int32)
        self.size = np.ones(n, dtype=np.int32)
        self.count = n
        self.largest = 1 if n else 0

    def find(self, x):
        x = int(x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = int(self.parent[x])
        return x

    def join(self, x, y):
        a, b = self.find(x), self.find(y)
        if a == b:
            return False
        if self.size[a] < self.size[b]:
            a, b = b, a
        self.parent[b] = a
        self.size[a] += self.size[b]
        self.count -= 1
        self.largest = max(self.largest, int(self.size[a]))
        return True


def asint(value):
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return None


def choose_fovs(metadata_path, samples):
    groups = {s[0]: Counter() for s in samples}
    with metadata_path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            fov = asint(row.get("fov"))
            if fov is None:
                continue
            for geo, _, start, end in samples:
                if start <= fov <= end:
                    groups[geo][fov] += 1
                    break
    selected = {}
    for geo, _, _, _ in samples:
        counts = groups[geo]
        if not counts:
            raise ValueError(f"No metadata FOVs found for {geo}")
        target = median(counts.values())
        # Representative by median cell count (not yet histological quality).
        selected[geo] = min(counts, key=lambda f: (abs(counts[f]-target), f))
    return selected, groups


def load_metadata(metadata_path, selected):
    per_fov = defaultdict(dict)
    with metadata_path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        columns = set(reader.fieldnames or [])
        required = {"fov", "cell", "CenterX_global_px", "CenterY_global_px"}
        if not required.issubset(columns):
            raise ValueError(f"Missing metadata columns: {required - columns}")
        for row in reader:
            fov = asint(row.get("fov"))
            if fov not in selected:
                continue
            key = (fov, row["cell"])
            try:
                x = float(row["CenterX_global_px"])
                y = float(row["CenterY_global_px"])
                area = float(row.get("Area") or "nan")
            except ValueError:
                continue
            if not all(map(math.isfinite, (x, y))):
                continue
            per_fov[fov][key] = (x, y, area)
    return per_fov


def load_expression(matrix_path, selected, meta, cutoff):
    wanted_genes = sorted({g for pair in PAIRS for g in pair[:2]})
    per_fov = defaultdict(list)
    unmatched = Counter()
    with gzip.open(matrix_path, "rt", newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)
        cols = {g: header.index(g) for g in wanted_genes if g in header}
        if len(cols) != len(wanted_genes):
            raise ValueError("Missing genes: " + str(set(wanted_genes) - set(cols)))
        i_fov, i_cell = header.index("fov"), header.index("cell")
        base = 3
        for row in reader:
            fov = asint(row[i_fov])
            if fov not in selected:
                continue
            key = (fov, row[i_cell])
            if key not in meta[fov]:
                unmatched[fov] += 1
                continue
            # Per-cell UMI QC. With 6519 genes and ~1k cells per FOV,
            # this streaming step is safe for RAM, but reads full input.
            total = sum(float(value or 0) for value in row[base:])
            if total < cutoff:
                continue
            vals = [float(row[cols[g]] or 0) for g in wanted_genes]
            per_fov[fov].append((key, vals, total))
    return per_fov, wanted_genes, unmatched


def candidate_edges(coords, values, genes, radius, threshold, max_pairs):
    n = len(coords)
    if n < 2:
        return np.empty((0, 2), dtype=np.int32), np.zeros(0, dtype=np.int16), 0
    pairs = cKDTree(coords).query_pairs(radius, output_type="ndarray")
    if len(pairs) > max_pairs:
        raise RuntimeError(f"Too many spatial pairs ({len(pairs):,}); choose a smaller FOV/radius.")
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    uv = np.zeros(len(pairs), dtype=np.int16)
    vu = np.zeros(len(pairs), dtype=np.int16)
    a, b = pairs[:, 0], pairs[:, 1]
    for ligand, receptor, _ in PAIRS:
        li, ri = gene_to_idx[ligand], gene_to_idx[receptor]
        uv += (values[a, li] >= threshold) & (values[b, ri] >= threshold)
        vu += (values[b, li] >= threshold) & (values[a, ri] >= threshold)
    strength = np.maximum(uv, vu)
    return pairs, strength, len(pairs)


def filtration(rows, strengths, n, max_k, min_k):
    """
    For k = max_k,..,min_k add spatial edges supported by >=k LR pairs.
    Vertices are present at birth 0; edge arrival at weight=1/k.
    Return snapshots and positive-length merge death values. Each remaining
    connected component has one infinite 0D interval at final threshold.
    """
    dsu = DSU(n)
    order = np.argsort(-strengths, kind="stable")
    counts = Counter()
    deaths = []
    last = 0
    snapshots = []
    for k in range(max_k, min_k - 1, -1):
        while last < len(order) and strengths[order[last]] >= k:
            j = order[last]
            if dsu.join(rows[j, 0], rows[j, 1]):
                deaths.append(1.0 / int(strengths[j]))
            last += 1
        snapshots.append({
            "min_supported_lr_pairs": k,
            "filtration_weight_max": 1.0 / k,
            "active_edges": last,
            "n_components": int(dsu.count),
            "largest_component_cells": int(dsu.largest),
            "largest_component_fraction": round(dsu.largest / max(n, 1), 6),
            "fraction_vertices_in_largest": round(dsu.largest / max(n, 1), 6),
        })
    return snapshots, deaths, int(dsu.count)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/pilot_6k")
    p.add_argument("--out", type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/outputs/h0_lr_pilot")
    p.add_argument("--radius-diameters", type=float, default=2.5,
        help="Spatial radius in units of median approximate cell diameter; default 2.5")
    p.add_argument("--min-count", type=float, default=1.0,
        help="Per-gene detected count needed to support a ligand/receptor; default 1")
    p.add_argument("--min-total-count", type=float, default=50.0,
        help="Per-cell UMI QC; default 50")
    p.add_argument("--max-cells", type=int, default=4000)
    p.add_argument("--max-pairs", type=int, default=1500000)
    p.add_argument("--min-k", type=int, default=1,
        help="Weakest LR-pair support considered in H0 filtration; default 1")
    args = p.parse_args()
    root, out = args.root.expanduser(), args.out.expanduser()
    if args.radius_diameters <= 0 or args.min_k < 1:
        p.error("radius-diameters must be >0 and min-k >= 1")
    out.mkdir(parents=True, exist_ok=True)

    summaries, curves, persistence = [], [], []
    for prefix, samples in PREFIXES.items():
        metadata = root/"metadata"/(prefix+"_metadata_file.csv")
        matrix = root/"counts"/(prefix+"_exprMat_file.csv.gz")
        if not metadata.exists() or not matrix.exists():
            p.error("Missing source file: " + str(metadata if not metadata.exists() else matrix))
        chosen, counts = choose_fovs(metadata, samples)
        wanted = set(chosen.values())
        meta = load_metadata(metadata, wanted)
        expressions, genes, unmatched = load_expression(
            matrix, wanted, meta, args.min_total_count)
        for geo, condition, _, _ in samples:
            fov = chosen[geo]
            records = expressions[fov]
            n = len(records)
            if n < 2 or n > args.max_cells:
                print(f"SKIPPING {geo}: QC cells {n}, allowed 2..{args.max_cells}")
                continue
            coords = np.array([meta[fov][rec[0]][:2] for rec in records], dtype=float)
            areas = [meta[fov][rec[0]][2] for rec in records]
            diameters = [2*math.sqrt(a / math.pi) for a in areas
                         if math.isfinite(a) and a > 0]
            if not diameters:
                raise RuntimeError(f"{geo}: no valid segmentation areas")
            med_diameter = median(diameters)
            radius_px = med_diameter * args.radius_diameters
            values = np.asarray([rec[1] for rec in records], dtype=np.float32)
            pair_indices, support, spatial_pairs = candidate_edges(
                coords, values, genes, radius_px, args.min_count, args.max_pairs)
            snapshots, deaths, infinite = filtration(
                pair_indices, support, n, len(PAIRS), args.min_k)
            geom_dsu = DSU(n)
            for a, b in pair_indices:
                geom_dsu.join(a, b)
            for point in snapshots:
                curves.append({"geo":geo,"condition":condition,"fov":fov, **point})
            for death in deaths:
                persistence.append({"geo":geo,"condition":condition,"fov":fov,
                                    "birth":0.0,"death":death})
            for _ in range(infinite):
                persistence.append({"geo":geo,"condition":condition,"fov":fov,
                                    "birth":0.0,"death":"inf"})
            summary = {
                "geo": geo, "condition": condition, "fov": fov,
                "fov_original_cell_count": counts[geo][fov],
                "n_cells_after_qc_and_join": n,
                "n_unmatched_expression_rows": unmatched[fov],
                "median_cell_diameter_px": round(med_diameter,3),
                "radius_px": round(radius_px,3),
                "radius_in_cell_diameters": args.radius_diameters,
                "n_geometric_neighbor_pairs": spatial_pairs,
                "geometric_components": geom_dsu.count,
                "n_edges_any_lr": int((support >= 1).sum()),
                "n_edges_two_or_more_lr": int((support >= 2).sum()),
                "n_infinite_h0_intervals_at_min_k": infinite,
                "panel_lr_pairs": len(PAIRS),
            }
            summaries.append(summary)
            print(f"{geo} {condition}, FOV {fov}: {n:,} QC cells; "
                  f"{spatial_pairs:,} spatial pairs; "
                  f"{summary['n_edges_any_lr']:,} LR-supported edges; "
                  f"geometric components {geom_dsu.count}; "
                  f"H0 infinite intervals {infinite}", flush=True)

    def write_csv(name, rows):
        if not rows:
            return
        with (out/name).open("w", newline="") as f:
            writer=csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    write_csv("sample_summary.csv", summaries)
    write_csv("component_filtration.csv", curves)
    write_csv("h0_intervals.csv", persistence)
    (out/"run_config.json").write_text(json.dumps({
        "parameters": vars(args) | {"root":str(root),"out":str(out)},
        "candidate_pairs":PAIRS,
        "interpretation":"Exploratory binary expression + spatial proximity H0. "
           "Neither COMMOT nor sheaf/receiver-response compatibility is calculated. "
           "Each FOV is a non-random biological region, not an independent donor. "
           "Weights 1/N are not physical concentration or calibrated potency."
    },indent=2))
    print(f"Results written to {out}. Do not infer disease effects from four FOVs.")

if __name__ == "__main__":
    main()
