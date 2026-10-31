#!/usr/bin/env python3
"""Full H0-filtration spatial-shuffle diagnostic for GSE312415 CosMx 6K.

Run this script next to h0_all_fovs.py (the earlier working analysis).
It uses the SAME LIANA resource, metadata joins, QC, spatial-radius rule,
ligand-complex support metric, and specimen/FOV mappings.

For each FOV, shuffle complete molecular profiles among fixed cell positions.
Compute H0 component count and giant-component fraction at EVERY observed or
null-positive ligand-support level. Saves a checkpoint per FOV so interrupted
runs can resume. No additional sequencing files are downloaded.

Scientific limits:
- This measures candidate communication POTENTIAL, not flux or sheaf compatibility.
- FOVs and cells are nested observations, NOT independent patients.
- Spatially shuffling cell profiles also shuffles cell types and any spatial
  coexpression patterns; it is not a composition-controlled causal null.
"""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from statistics import median

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import h0_all_fovs as base

SHARED_THRESHOLDS = [1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 75, 100, 150, 200]
DEFAULT_RESOURCE_SHA256 = "65e46a428818c0600de66c1574c94ff2f5dc80859998321c03a956d330dd67fb"


def evaluate_filtration(n, edges, strength, thresholds):
    """Exact union-find filtration, recording H0 at every requested k."""
    good = np.flatnonzero(strength > 0)
    order = good[np.argsort(-strength[good], kind="stable")]
    dsu = base.DSU(n)
    components = np.empty(len(thresholds), dtype=np.int32)
    giants = np.empty(len(thresholds), dtype=np.float32)
    edge_counts = np.empty(len(thresholds), dtype=np.int32)
    ptr = 0
    for i, k in enumerate(thresholds):
        while ptr < len(order) and strength[order[ptr]] >= k:
            row = edges[order[ptr]]
            dsu.union(int(row[0]), int(row[1]))
            ptr += 1
        components[i] = dsu.components
        giants[i] = dsu.largest / n if n else 0
        edge_counts[i] = ptr
    return components, giants, edge_counts


def edge_support(u, v, X, Y, batch=6000):
    """Maximum of u->v and v->u supported-pathway counts."""
    support = np.empty(len(u), dtype=np.int32)
    for i in range(0, len(u), batch):
        j = min(i + batch, len(u))
        forward = np.asarray(X[u[i:j]].multiply(Y[v[i:j]]).sum(axis=1)).ravel()
        reverse = np.asarray(X[v[i:j]].multiply(Y[u[i:j]]).sum(axis=1)).ravel()
        support[i:j] = np.maximum(forward, reverse).astype(np.int32)
    return support


def process_fov(cell_records, meta_fov, res, radius_multiplier,
                metric, permutations, seed, max_pairs, max_cells):
    n = len(cell_records)
    if n < 2 or n > max_cells:
        return None
    points = np.array([meta_fov[cell][:2] for cell, _ in cell_records], dtype=np.float64)
    areas = np.array([meta_fov[cell][2] for cell, _ in cell_records], dtype=np.float64)
    diameter_px = float(np.median(2.0 * np.sqrt(areas / np.pi)))
    radius_px = diameter_px * radius_multiplier
    pairs = cKDTree(points).query_pairs(radius_px, output_type="ndarray")
    if len(pairs) > max_pairs:
        raise RuntimeError(f"Too many geometric neighbour pairs: {len(pairs):,}")
    pairs = pairs.astype(np.int32).reshape((-1, 2))
    if len(pairs):
        pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]
    u, v = pairs[:, 0], pairs[:, 1]

    binary = np.stack([sig for _, sig in cell_records])
    L = base.make_complex_matrix(binary, res["needed_genes"], res["ligands"])
    R = base.make_complex_matrix(binary, res["needed_genes"], res["receptors"])
    if metric == "ligand":
        # Treat a ligand complex as supported if any compatible receptor
        # complex is detected. Each ligand counts at most once per edge direction.
        matched_receptors = (R @ res["adjacency"].T).tocsr()
        matched_receptors.data[:] = 1
        matched_receptors.eliminate_zeros()
        first, second = L, matched_receptors
    else:
        first, second = (L @ res["adjacency"]).tocsr(), R

    observed_support = edge_support(u, v, first, second)
    rng = np.random.default_rng(seed)
    null_supports = []
    for _ in range(permutations):
        perm = rng.permutation(n)
        null_supports.append(edge_support(perm[u], perm[v], first, second))

    levels = set(SHARED_THRESHOLDS)
    levels.update(int(x) for x in np.unique(observed_support) if x > 0)
    for strengths in null_supports:
        levels.update(int(x) for x in np.unique(strengths) if x > 0)
    thresholds = np.array(sorted(levels, reverse=True), dtype=np.int32)
    obs_components, obs_giant, obs_edges = evaluate_filtration(
        n, pairs, observed_support, thresholds)
    null_components = np.empty((permutations, len(thresholds)), dtype=np.int32)
    null_giant = np.empty((permutations, len(thresholds)), dtype=np.float32)
    null_edges = np.empty((permutations, len(thresholds)), dtype=np.int32)
    for j, strengths in enumerate(null_supports):
        c, g, e = evaluate_filtration(n, pairs, strengths, thresholds)
        null_components[j] = c
        null_giant[j] = g
        null_edges[j] = e

    return {
        "n_cells": np.array(n, dtype=np.int32),
        "spatial_pairs": np.array(len(pairs), dtype=np.int32),
        "median_diameter_px": np.array(diameter_px),
        "radius_px": np.array(radius_px),
        "thresholds": thresholds,
        "observed_components": obs_components,
        "observed_giant": obs_giant,
        "observed_edges": obs_edges,
        "null_components": null_components,
        "null_giant": null_giant,
        "null_edges": null_edges,
    }


def write_csv(path, records, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=fields)
        w.writeheader()
        w.writerows(records)


def summarise_checkpoints(checkpoint_dir, out, permutations):
    records = []
    for path in sorted(checkpoint_dir.glob("*.npz")):
        name = path.stem
        geo, fovstr = name.rsplit("_fov", 1)
        fov = int(fovstr)
        condition = next((cond for group in base.SPECIMENS.values()
                          for id_, cond, _, _ in group if id_ == geo), "unknown")
        with np.load(path, allow_pickle=False) as data:
            thresholds = data["thresholds"]
            n = int(data["n_cells"])
            nc = data["null_components"]
            ng = data["null_giant"]
            ne = data["null_edges"]
            if nc.shape[0] != permutations:
                raise RuntimeError(
                    f"{name} saved with {nc.shape[0]} null permutations, "
                    f"requested {permutations}; use a new --out directory")
            for j, k in enumerate(thresholds):
                c = int(data["observed_components"][j])
                g = float(data["observed_giant"][j])
                e = int(data["observed_edges"][j])
                cmean = float(nc[:, j].mean())
                gmean = float(ng[:, j].mean())
                records.append({
                    "geo": geo, "condition": condition, "fov": fov,
                    "n_cells": n, "min_support": int(k),
                    "observed_components": c,
                    "null_mean_components": round(cmean, 5),
                    "delta_components_per_cell": round((c-cmean)/n, 7),
                    "observed_largest_fraction": round(g, 7),
                    "null_mean_largest_fraction": round(gmean, 7),
                    "delta_largest_fraction": round(g-gmean, 7),
                    "null_sd_largest_fraction": round(
                        float(ng[:, j].std(ddof=1)), 7),
                    "observed_edges": e,
                    "null_mean_edges": round(float(ne[:, j].mean()), 3),
                    "null_sd_components": round(float(nc[:, j].std(ddof=1)), 4)
                })
    if not records:
        print("No checkpoints available yet.", flush=True)
        return
    write_csv(out/"full_filtration_vs_null.csv", records, list(records[0]))
    grouped = {}
    for row in records:
        if row["min_support"] not in SHARED_THRESHOLDS:
            continue
        grouped.setdefault((row["geo"], row["condition"], row["min_support"]), []).append(row)
    summary = []
    for (geo, condition, k), subset in sorted(grouped.items()):
        summary.append({
            "geo": geo, "condition": condition, "min_support": k,
            "n_fovs": len(subset),
            "median_observed_largest_fraction": round(median(
                x["observed_largest_fraction"] for x in subset), 6),
            "median_null_largest_fraction": round(median(
                x["null_mean_largest_fraction"] for x in subset), 6),
            "median_fov_delta_largest_fraction": round(median(
                x["delta_largest_fraction"] for x in subset), 6),
            "median_fov_delta_components_per_cell": round(median(
                x["delta_components_per_cell"] for x in subset), 6),
        })
    write_csv(out/"specimen_threshold_summary.csv", summary, list(summary[0]))
    print(f"Summarised {len(list(checkpoint_dir.glob('*.npz')))} FOVs into "
          f"{len(records)} observed-null threshold comparisons.", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/pilot_6k")
    ap.add_argument("--out", type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/outputs/h0_full_null")
    ap.add_argument("--catalog", type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/metadata/liana_omni_resource.csv")
    ap.add_argument("--resource", default="consensus")
    ap.add_argument("--metric", choices=["ligand", "pair"], default="ligand")
    ap.add_argument("--radius-diameters", type=float, default=2.5)
    ap.add_argument("--detection-count", type=float, default=1.0)
    ap.add_argument("--min-total-count", type=float, default=50.0)
    ap.add_argument("--permutations", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20261031)
    ap.add_argument("--max-cells-per-fov", type=int, default=6000)
    ap.add_argument("--max-neighbor-pairs", type=int, default=1500000)
    ap.add_argument("--chunk-size", type=int, default=384)
    ap.add_argument("--max-fovs", type=int, default=0,
                    help="For initial smoke test: process at most this many FOVs (0 = all)")
    args = ap.parse_args()
    if args.permutations < 2 or args.radius_diameters <= 0:
        ap.error("permutations must be >=2 and radius-diameters >0")
    out = args.out.expanduser()
    root = args.root.expanduser()
    catalog = args.catalog.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = out/"checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)
    cat_sha = base.get_catalog(catalog)
    if cat_sha != DEFAULT_RESOURCE_SHA256:
        raise RuntimeError(f"LIANA catalogue SHA256 changed: {cat_sha}. "
                           "Confirm the source version before reproducing this experiment.")
    config = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
              if k not in ("max_fovs",)}
    config.update({"catalog_sha256":cat_sha,
                   "analysis":"H0 on spatial signalling potential, not COMMOT or sheaf.",
                   "null":"Random permutation of entire cell expression profiles within each FOV.",
                   "independence":"FOVs nested within four GEO specimens; not independent patients."})
    config_path = out/"run_config.json"
    if config_path.exists():
        old = json.loads(config_path.read_text())
        if old != config:
            raise RuntimeError(f"Existing {out} has different settings. "
                               "Choose a new --out directory to avoid mixing results.")
    else:
        config_path.write_text(json.dumps(config, indent=2))
    rng = np.random.default_rng(args.seed)
    visited = 0
    completed = 0

    for prefix, specimens in base.SPECIMENS.items():
        expr_path = root/"counts"/(prefix+"_exprMat_file.csv.gz")
        metadata_path = root/"metadata"/(prefix+"_metadata_file.csv")
        if not expr_path.exists() or not metadata_path.exists():
            raise FileNotFoundError(f"Missing {expr_path} or {metadata_path}")
        meta, fov_to_geo = base.load_metadata(metadata_path, specimens)
        pending = [fov for fov in sorted(meta) if fov_to_geo.get(fov) and
                   not (checkpoint_dir / f"{fov_to_geo[fov]}_fov{fov}.npz").exists()]
        if not pending:
            print(f"{prefix}: already checkpointed, skipping expression reload.", flush=True)
            continue
        if args.max_fovs and visited >= args.max_fovs:
            break
        header = pd.read_csv(expr_path, compression="gzip", nrows=0)
        genes = set(map(str.upper, header.columns[3:]))
        resource = base.catalog_complexes(catalog, genes, args.resource)
        print(f"{prefix}: {len(pending)} FOVs still pending; "
              f"{len(resource['pairs'])} LR interactions.", flush=True)
        selected_meta = {fov: meta[fov] for fov in pending}
        records, matched, passing = base.load_counts(
            expr_path, selected_meta, resource["needed_genes"],
            args.detection_count, args.min_total_count, args.chunk_size)
        for fov in pending:
            if args.max_fovs and visited >= args.max_fovs:
                break
            visited += 1
            geo = fov_to_geo[fov]
            cells = records.get(fov, [])
            # Independent seed per specimen/FOV is stable across resumed runs.
            fov_seed = int.from_bytes(
                hashlib.sha256(f"{args.seed}:{geo}:{fov}".encode()).digest()[:8],
                byteorder="big")
            result = process_fov(
                cells, meta[fov], resource, args.radius_diameters,
                args.metric, args.permutations, fov_seed,
                args.max_neighbor_pairs, args.max_cells_per_fov)
            if result is None:
                print(f"{geo} FOV {fov}: skipped, {len(cells)} QC cells.", flush=True)
                continue
            dest = checkpoint_dir/f"{geo}_fov{fov}.npz"
            tmp = dest.with_suffix(".tmp")
            with tmp.open("wb") as handle:
                np.savez_compressed(handle, **result)
            os.replace(tmp, dest)
            completed += 1
            idx = np.where(result["thresholds"] == 20)[0]
            effect = (float(result["observed_giant"][idx[0]]) -
                      float(result["null_giant"][:, idx[0]].mean())) if len(idx) else float("nan")
            print(f"{geo} FOV {fov}: n={int(result['n_cells']):,}, "
                  f"pairs={int(result['spatial_pairs']):,}; "
                  f"observed-minus-null giant at k=20: {effect:+.3f}", flush=True)

    summarise_checkpoints(checkpoint_dir, out, args.permutations)
    print(f"New FOVs completed: {completed}. Outputs in {out}", flush=True)


if __name__ == "__main__":
    main()
