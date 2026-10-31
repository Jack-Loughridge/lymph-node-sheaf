#!/usr/bin/env python3
"""H0 full-filtration nulls: unrestricted versus broad-cell-type-preserving shuffles.

GSE312415 CosMx 6K, four specimens. Reuses h0_all_fovs.py for
FOV/specimen separation, ligand/receptor resource, coordinates, and RNA QC.

IMPORTANT
- Spatial support is candidate signalling potential, NOT sheaf compatibility.
- An automated marker-derived class is a PROXY, not a validated cell annotation.
- Cell/FOV observations are nested in specimens, not independent patients.
- Full profile permutations preserve the set of expression profiles per FOV;
  within-class permutations additionally preserve the broad class spatial map.
- No raw or individual-cell data are written to GitHub.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
from statistics import median
import sys

import numpy as np
import pandas as pd
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

import h0_all_fovs as base

EXPECTED_CATALOG_SHA = "65e46a428818c0600de66c1574c94ff2f5dc80859998321c03a956d330dd67fb"
MARKERS = {
    "epithelial": ("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CLDN4", "MUC2"),
    "t_cell": ("CD3D", "CD3E", "TRAC", "CD2", "CD247"),
    "b_cell": ("MS4A1", "CD79A", "CD79B", "CD19", "CD22"),
    "plasma_cell": ("MZB1", "JCHAIN", "XBP1", "SDC1", "IGKC"),
    "myeloid": ("LST1", "CD68", "LYZ", "FCGR3A", "CSF1R"),
    "fibroblast": ("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"),
    "endothelial": ("PECAM1", "VWF", "KDR", "EMCN", "CLDN5"),
    "mast_cell": ("TPSAB1", "TPSB2", "CPA3", "KIT"),
    "pericyte": ("RGS5", "PDGFRB", "CSPG4", "MCAM"),
    "smooth_muscle": ("ACTA2", "TAGLN", "MYH11"),
}
LABEL_COLUMNS = [
    "cell_type", "celltype", "cell_type_annotation", "predicted_cell_type",
    "major_cell_type", "broad_cell_type", "cell_type_major",
    "annotation", "celltype_l1", "celltype_l2"
]
SHARED_THRESHOLDS = [1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 75, 100, 150, 200]


def find_label_column(path, user_col):
    with path.open(newline="", encoding="utf-8-sig") as f:
        cols = csv.DictReader(f).fieldnames or []
    if user_col:
        if user_col not in cols:
            raise ValueError(f"Requested --label-column {user_col!r} missing; columns: {cols}")
        return user_col
    lower = {c.lower(): c for c in cols}
    for candidate in LABEL_COLUMNS:
        if candidate in lower:
            return lower[candidate]
    return None


def read_labels(path, column, wanted_fovs):
    if column is None:
        return {}
    labels = {}
    with path.open(newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            fov = base.as_fov(row.get("fov"))
            if fov not in wanted_fovs:
                continue
            label = str(row.get(column, "") or "").strip()
            if label.lower() in ("", "na", "nan", "none", "unknown", "unassigned"):
                label = "unclassified"
            labels[(fov, str(row["cell"]))] = label
    return labels


def assign_marker_groups(binary, gene_names):
    """Broad exploratory proxy using presence/absence of canonical markers.

    Require >=2 markers of one lineage, a >=0.25 fraction of panel markers,
    and a >=0.10 margin over the next lineage. Others are unclassified.
    """
    g_index = {g: j for j, g in enumerate(gene_names)}
    names, score, detected = [], [], []
    coverage = {}
    for name, markers in MARKERS.items():
        covered = [g for g in markers if g in g_index]
        coverage[name] = covered
        if len(covered) < 2:
            continue
        indices = [g_index[g] for g in covered]
        n_found = binary[:, indices].sum(axis=1)
        names.append(name)
        score.append(n_found / len(indices))
        detected.append(n_found)
    if not names:
        return np.array(["unclassified"] * len(binary)), coverage
    scores = np.column_stack(score)
    counts = np.column_stack(detected)
    best = np.argmax(scores, axis=1)
    best_score = scores[np.arange(len(binary)), best]
    best_count = counts[np.arange(len(binary)), best]
    if scores.shape[1] > 1:
        second_score = np.partition(scores, -2, axis=1)[:, -2]
    else:
        second_score = np.zeros(len(binary))
    valid = ((best_count >= 2) & (best_score >= 0.25)
             & (best_score - second_score >= 0.10))
    classes = np.array(["unclassified"] * len(binary), dtype=object)
    names = np.array(names, dtype=object)
    classes[valid] = names[best[valid]]
    return classes.astype(str), coverage


def get_cell_groups(records, required_genes, label_map, fov,
                    label_col, mode):
    if mode != "markers" and label_col is not None:
        classes = np.asarray([
            label_map.get((fov, cell), "unclassified")
            for cell, _ in records], dtype=str)
        assigned = np.mean(classes != "unclassified")
        if mode == "metadata" or (assigned >= 0.8 and len(set(classes)) >= 2):
            return classes, "metadata:" + label_col, {}
        print(f"  FOV {fov}: metadata field {label_col} only "
              f"{assigned:.1%} assigned; using marker proxies.", flush=True)
    if mode == "metadata":
        raise RuntimeError("Cell-type metadata requested but not available")
    signatures = np.stack([b for _, b in records])
    classes, covered = assign_marker_groups(signatures, required_genes)
    return classes, "marker_proxies", covered


def edge_support(u, v, first, second, batch=5000):
    """Count distinct detected ligand complexes, or LR complexes if chosen."""
    result = np.empty(len(u), dtype=np.int32)
    for start in range(0, len(u), batch):
        stop = min(start + batch, len(u))
        fwd = np.asarray(
            first[u[start:stop]].multiply(second[v[start:stop]]).sum(axis=1)
        ).ravel()
        rev = np.asarray(
            first[v[start:stop]].multiply(second[u[start:stop]]).sum(axis=1)
        ).ravel()
        result[start:stop] = np.maximum(fwd, rev).astype(np.int32)
    return result


def forest_curves(n, pairs, strength, max_k):
    """Exact H0 curve from maximum-support spanning forest.

    Weighted minimum spanning forest on cost=1/support retains the SAME
    connected components at every k as thresholding ALL edges at support>=k.
    The forest contains <= n-1 edges, greatly accelerating 100s of nulls.
    """
    positive = strength > 0
    p = pairs[positive]
    w = strength[positive]
    counts = np.bincount(w, minlength=max_k + 1)
    active_at_k = np.cumsum(counts[::-1])[::-1]

    forest = base.DSU(n)
    comps = np.empty(max_k + 1, dtype=np.int32)
    giant = np.empty(max_k + 1, dtype=np.float32)
    comps[0], giant[0] = n, 1.0 / max(n, 1)
    if len(p):
        cost = 1.0 / w.astype(np.float64)
        symmetric = csr_matrix(
            (np.concatenate((cost, cost)),
             (np.concatenate((p[:, 0], p[:, 1])),
              np.concatenate((p[:, 1], p[:, 0])))),
            shape=(n, n))
        mst = minimum_spanning_tree(symmetric).tocoo()
        mstrength = np.rint(1.0 / mst.data).astype(np.int32)
        order = np.argsort(-mstrength, kind="stable")
        mr, mc = mst.row[order], mst.col[order]
        mw = mstrength[order]
    else:
        mr = mc = mw = np.empty(0, dtype=np.int32)
    idx = 0
    for k in range(max_k, 0, -1):
        while idx < len(mw) and mw[idx] >= k:
            forest.union(int(mr[idx]), int(mc[idx]))
            idx += 1
        comps[k] = forest.components
        giant[k] = forest.largest / max(n, 1)
    return comps, giant, active_at_k


def permutation(n, kind, classes, rng):
    if kind == "unrestricted":
        return rng.permutation(n)
    perm = np.arange(n)
    for label in np.unique(classes):
        indices = np.flatnonzero(classes == label)
        if len(indices) > 1:
            perm[indices] = rng.permutation(indices)
    return perm


def analysis_fov(records, meta_fov, resource, metric, radius_mult,
                 permutations, seed, max_cells, max_pairs,
                 type_labels, label_source):
    n = len(records)
    if n < 2 or n > max_cells:
        return None
    points = np.asarray([meta_fov[cell][:2] for cell, _ in records],
                        dtype=np.float64)
    areas = np.asarray([meta_fov[cell][2] for cell, _ in records],
                       dtype=np.float64)
    diameter = float(np.median(2 * np.sqrt(areas / np.pi)))
    radius = diameter * radius_mult
    pairs = cKDTree(points).query_pairs(radius, output_type="ndarray")
    if len(pairs) > max_pairs:
        raise RuntimeError(f"Too many neighbour pairs: {len(pairs)}")
    pairs = pairs.astype(np.int32).reshape(-1, 2)
    u, v = pairs[:, 0], pairs[:, 1]
    matrix = np.stack([b for _, b in records])
    L = base.make_complex_matrix(
        matrix, resource["needed_genes"], resource["ligands"])
    R = base.make_complex_matrix(
        matrix, resource["needed_genes"], resource["receptors"])
    if metric == "ligand":
        Q = (R @ resource["adjacency"].T).tocsr()
        Q.data[:] = 1
        Q.eliminate_zeros()
        first, second = L, Q
    else:
        first, second = (L @ resource["adjacency"]).tocsr(), R

    observed = edge_support(u, v, first, second)
    rng = np.random.default_rng(seed)
    nulls = {"unrestricted": [], "within_type": []}
    moved = {"unrestricted": [], "within_type": []}
    max_support = max(200, int(observed.max(initial=0)))
    for kind in nulls:
        for j in range(permutations):
            perm = permutation(n, kind, type_labels, rng)
            moved[kind].append(float(np.mean(perm != np.arange(n))))
            arr = edge_support(perm[u], perm[v], first, second)
            nulls[kind].append(arr)
            max_support = max(max_support, int(arr.max(initial=0)))

    oc, og, oe = forest_curves(n, pairs, observed, max_support)
    output = {
        "n_cells": np.array(n, dtype=np.int32),
        "n_spatial_pairs": np.array(len(pairs), dtype=np.int32),
        "max_support": np.array(max_support, dtype=np.int32),
        "median_diameter_px": np.array(diameter),
        "radius_px": np.array(radius),
        "observed_components": oc,
        "observed_giant": og,
        "observed_edges": oe,
        "type_assigned_fraction": np.array(
            float(np.mean(type_labels != "unclassified"))),
        "type_source": np.array(label_source),
        "type_counts_json": np.array(json.dumps(dict(Counter(type_labels)))),
    }
    for kind in nulls:
        c_rows, g_rows, e_rows = [], [], []
        for weights in nulls[kind]:
            c, g, e = forest_curves(n, pairs, weights, max_support)
            c_rows.append(c)
            g_rows.append(g)
            e_rows.append(e)
        output[kind+"_components"] = np.stack(c_rows)
        output[kind+"_giant"] = np.stack(g_rows)
        output[kind+"_edges"] = np.stack(e_rows)
        output[kind+"_moved"] = np.asarray(moved[kind])
    return output


def write_csv(path, data):
    if not data:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0]), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)


def summarize_checkpoints(directory, out, n_permutations):
    rows = []
    fov_summary = []
    for path in sorted(directory.glob("GSM*_fov*.npz")):
        geo, fovstr = path.stem.rsplit("_fov", 1)
        fov = int(fovstr)
        condition = next((cond for group in base.SPECIMENS.values()
                          for g, cond, lo, hi in group if g == geo),
                         "unknown")
        with np.load(path, allow_pickle=False) as d:
            n = int(d["n_cells"])
            maxk = int(d["max_support"])
            groups = json.loads(str(d["type_counts_json"]))
            fov_summary.append({
                "geo": geo, "condition": condition, "fov": fov,
                "n_cells": n, "n_spatial_pairs": int(d["n_spatial_pairs"]),
                "type_source": str(d["type_source"]),
                "type_assigned_fraction": float(d["type_assigned_fraction"]),
                "type_counts_json": json.dumps(groups),
                "within_type_fraction_moved": float(d["within_type_moved"].mean()),
                "unrestricted_fraction_moved": float(d["unrestricted_moved"].mean()),
            })
            for k in range(1, maxk + 1):
                g = float(d["observed_giant"][k])
                c = int(d["observed_components"][k])
                row = {
                    "geo":geo,"condition":condition,"fov":fov,
                    "n_cells":n,"min_support":k,
                    "observed_largest_fraction":g,
                    "observed_components_per_cell":c / n,
                    "observed_edges":int(d["observed_edges"][k])
                }
                for kind in ("unrestricted", "within_type"):
                    ng = d[kind+"_giant"][:, k]
                    nc = d[kind+"_components"][:, k] / n
                    ne = d[kind+"_edges"][:, k]
                    if len(ng) != n_permutations:
                        raise RuntimeError(
                            f"{path.name} has {len(ng)} {kind} permutations, "
                            f"requested {n_permutations}")
                    row[kind+"_null_mean_largest"] = float(ng.mean())
                    row[kind+"_delta_largest"] = g - float(ng.mean())
                    row[kind+"_null_sd_largest"] = float(ng.std(ddof=1))
                    row[kind+"_null_mean_components_per_cell"] = float(nc.mean())
                    row[kind+"_delta_components_per_cell"] = c / n - float(nc.mean())
                    row[kind+"_null_mean_edges"] = float(ne.mean())
                    # Descriptive per-FOV permutation tail rates. No FOVs
                    # treated as independent participants.
                    row[kind+"_p_giant_high"] = (1+int(np.sum(ng>=g-1e-8))) / (n_permutations+1)
                    row[kind+"_p_giant_low"] = (1+int(np.sum(ng<=g+1e-8))) / (n_permutations+1)
                rows.append(row)
    if not rows:
        print("No completed checkpoints yet.", flush=True)
        return
    write_csv(out/"full_filtration_vs_two_nulls.csv", rows)
    write_csv(out/"fov_type_summary.csv", fov_summary)
    by_spec = defaultdict(list)
    for r in rows:
        if r["min_support"] in SHARED_THRESHOLDS:
            by_spec[(r["geo"], r["condition"], r["min_support"])].append(r)
    aggregated = []
    for (geo, cond, k), sub in sorted(by_spec.items()):
        aggregated.append({
            "geo":geo, "condition":cond, "min_support":k,
            "n_fovs":len(sub),
            "median_observed_largest":float(median(r["observed_largest_fraction"] for r in sub)),
            "median_unrestricted_delta":float(median(r["unrestricted_delta_largest"] for r in sub)),
            "median_within_type_delta":float(median(r["within_type_delta_largest"] for r in sub)),
            "median_within_type_components_delta":float(median(r["within_type_delta_components_per_cell"] for r in sub)),
        })
    write_csv(out/"specimen_threshold_summary.csv", aggregated)
    print(f"Summarised {len(fov_summary)} FOVs and {len(rows)} thresholds. "
          f"Outputs: {out}", flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/pilot_6k")
    p.add_argument("--catalog", type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/metadata/liana_omni_resource.csv")
    p.add_argument("--out", type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/outputs/h0_stratified_100")
    p.add_argument("--resource", choices=("consensus","cellchatdb","cellphonedb"),
                   default="consensus")
    p.add_argument("--metric", choices=("ligand","pair"), default="ligand")
    p.add_argument("--type-source", choices=("auto","metadata","markers"),
                   default="auto")
    p.add_argument("--label-column", default=None)
    p.add_argument("--radius-diameters", type=float, default=2.5)
    p.add_argument("--detection-count", type=float, default=1.0)
    p.add_argument("--min-total-count", type=float, default=50.0)
    p.add_argument("--permutations", type=int, default=100)
    p.add_argument("--seed", type=int, default=20261031)
    p.add_argument("--max-fovs", type=int, default=0,
                   help="Limit NEW FOVs processed in this invocation; 0=all")
    p.add_argument("--max-cells-per-fov", type=int, default=6000)
    p.add_argument("--max-neighbor-pairs", type=int, default=1500000)
    p.add_argument("--chunk-size", type=int, default=384)
    args = p.parse_args()
    if args.permutations < 2 or args.radius_diameters <= 0 or args.max_fovs < 0:
        p.error("Need >=2 permutations, radius >0 and max-fovs >=0")
    out, root, cat = args.out.expanduser(), args.root.expanduser(), args.catalog.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    checkpoints = out/"checkpoints"
    checkpoints.mkdir(exist_ok=True)
    sha = base.get_catalog(cat)
    if sha != EXPECTED_CATALOG_SHA:
        raise RuntimeError("LIANA resource changed; inspect before using: "+sha)
    config = {k: str(v) if isinstance(v,Path) else v for k,v in vars(args).items()
              if k != "max_fovs"}
    config.update({"catalog_sha256":sha,
                   "description":"Potential-signalling H0, full threshold range, "
                       "unrestricted and within-broad-cell-type permutations.",
                   "proxy_caveat":"Markers produce provisional broad classes, "
                       "not validated cell types."})
    cfg = out/"run_config.json"
    if cfg.exists() and json.loads(cfg.read_text()) != config:
        raise RuntimeError("Existing output directory contains different parameters. "
                           "Use a fresh --out when changing settings.")
    cfg.write_text(json.dumps(config,indent=2))

    done=0
    for prefix, specimens in base.SPECIMENS.items():
        expr_file=root/"counts"/(prefix+"_exprMat_file.csv.gz")
        meta_file=root/"metadata"/(prefix+"_metadata_file.csv")
        meta, fov_to_geo = base.load_metadata(meta_file, specimens)
        pending = [f for f in sorted(meta)
                   if not (checkpoints/f"{fov_to_geo[f]}_fov{f}.npz").exists()]
        if not pending:
            print(f"{prefix}: complete, skipping counts reload.", flush=True)
            continue
        if args.max_fovs and done >= args.max_fovs:
            break
        col = None if args.type_source == "markers" else find_label_column(
            meta_file,args.label_column)
        if args.type_source == "metadata" and col is None:
            raise RuntimeError(f"No cell-type annotation column in {meta_file}")
        labels_map = read_labels(meta_file,col,set(pending))
        panel = pd.read_csv(expr_file,compression="gzip",nrows=0)
        panel_genes=set(str(g).upper() for g in panel.columns[3:])
        resource=base.catalog_complexes(cat,panel_genes,args.resource)
        needed=sorted(set(resource["needed_genes"]) |
                      {g for group in MARKERS.values() for g in group
                       if g in panel_genes})
        resource["needed_genes"]=needed
        print(f"{prefix}: pending FOVs {len(pending)}; "
              f"{len(resource['pairs'])} interactions, "
              f"{len(needed)} relevant gene columns; labels={col or 'marker proxies'}",
              flush=True)
        records, matched, passing=base.load_counts(
            expr_file,{f:meta[f] for f in pending},needed,
            args.detection_count,args.min_total_count,args.chunk_size)
        for fov in pending:
            if args.max_fovs and done>=args.max_fovs:
                break
            geo=fov_to_geo[fov]
            cells=records.get(fov,[])
            if len(cells)<2 or len(cells)>args.max_cells_per_fov:
                print(f"Skipping {geo} FOV {fov}: {len(cells)} cells",flush=True)
                continue
            types, source, coverage=get_cell_groups(
                cells,needed,labels_map,fov,col,args.type_source)
            # Per-FOV reproducible independent random seeds.
            seed=int.from_bytes(hashlib.sha256(
                f"{args.seed}/{geo}/{fov}/stratified-v1".encode()).digest()[:8],"big")
            result=analysis_fov(
                cells,meta[fov],resource,args.metric,args.radius_diameters,
                args.permutations,seed,args.max_cells_per_fov,
                args.max_neighbor_pairs,types,source)
            dest=checkpoints/f"{geo}_fov{fov}.npz"
            tmp=dest.with_suffix(".part")
            with tmp.open("wb") as handle:
                np.savez_compressed(handle,**result)
            os.replace(tmp,dest)
            done+=1
            k=min(30,int(result["max_support"]))
            observed=float(result["observed_giant"][k])
            unrestricted=float(result["unrestricted_giant"][:,k].mean())
            stratified=float(result["within_type_giant"][:,k].mean())
            print(f"{geo} FOV {fov}: cells={len(cells)}, types={source}, "
                  f"labelled={float(result['type_assigned_fraction']):.1%}, "
                  f"moved(type)={result['within_type_moved'].mean():.1%}, "
                  f"k={k} delta(all)={observed-unrestricted:+.3f}, "
                  f"delta(type)={observed-stratified:+.3f}",flush=True)
    summarize_checkpoints(checkpoints,out,args.permutations)
    print("NOTE: compare specimens descriptively; only four independent specimens.",
          flush=True)


if __name__ == "__main__":
    main()
