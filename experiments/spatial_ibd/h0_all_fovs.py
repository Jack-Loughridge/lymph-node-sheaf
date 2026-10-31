#!/usr/bin/env python3
"""Whole-FOV spatial ligand-receptor H0 diagnostic for GEO GSE312415 CosMx 6K.

No sheaf, COMMOT or physical signal transport is estimated here.
Cell-centroid proximity + expression of curated ligand/receptor COMPLEXES
defines possible directed interactions. H0 is computed on the undirected
OR of the two directions, weighted by count of distinct supported ligands.

All FOVs are analysed, but FOVs and cells are NOT independent patients.
Two expression files each contain two separate specimens. Never connect FOVs
across specimens. Refer to run_config.json for all definitions.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from statistics import median
import sys
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

CATALOG_URL = ("https://raw.githubusercontent.com/saezlab/liana-py/"
               "main/src/liana/resource/omni_resource.csv")
SPECIMENS = {
    "Cos1400R355RRNA6K4u1716": [
        ("GSM9345416", "Rectum_HC", 32, 65),
        ("GSM9345417", "Rectum_UC_Noninflamed", 1, 31)],
    "Cos1380R387RRNA6K4u2716": [
        ("GSM9345418", "Rectum_UC_Inflamed", 18, 45),
        ("GSM9345419", "Rectum_UC_Inflamed", 1, 17)],
}
VALID_SYMBOL = re.compile(r"^[A-Za-z0-9_-]+$")


class DSU:
    def __init__(self, n):
        self.parent = np.arange(n, dtype=np.int32)
        self.sizes = np.ones(n, dtype=np.int32)
        self.components = n
        self.largest = int(n > 0)

    def find(self, i):
        i = int(i)
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = int(self.parent[i])
        return i

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return False
        if self.sizes[a] < self.sizes[b]:
            a, b = b, a
        self.parent[b] = a
        self.sizes[a] += self.sizes[b]
        self.components -= 1
        self.largest = max(self.largest, int(self.sizes[a]))
        return True


def get_catalog(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        print("Downloading LIANA consensus resource catalogue...", flush=True)
        req = Request(CATALOG_URL, headers={"User-Agent": "spatial-ibd-h0-audit/1.0"})
        tmp = path.with_name(path.name + ".part")
        try:
            with urlopen(req, timeout=90) as src, tmp.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            tmp.replace(path)
        except Exception:
            tmp.unlink(missing_ok=True)
            raise
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if len(raw) < 1000:
        raise ValueError("Ligand-receptor catalogue unexpectedly small")
    return digest


def catalog_complexes(path, genes, resource):
    # A complex such as ITGA1_ITGB1 is supported only if ALL components
    # are detectable. Do not split complexes into independent LR pairs.
    all_pairs = set()
    resource_rows = 0
    invalid = Counter()
    with path.open(newline="", encoding="utf-8-sig") as f:
        r = csv.DictReader(f)
        needed = {"source_genesymbol", "target_genesymbol", "resource"}
        if not needed.issubset(set(r.fieldnames or [])):
            raise ValueError(f"Wrong LIANA resource format: {r.fieldnames}")
        for row in r:
            if row["resource"].strip().lower() != resource:
                continue
            resource_rows += 1
            ligand = row["source_genesymbol"].strip().upper()
            receptor = row["target_genesymbol"].strip().upper()
            if not ligand or not receptor or not VALID_SYMBOL.fullmatch(ligand) or not VALID_SYMBOL.fullmatch(receptor):
                invalid["invalid_syntax"] += 1
                continue
            li = tuple(ligand.split("_"))
            ri = tuple(receptor.split("_"))
            if len(li) > 8 or len(ri) > 8 or "" in li or "" in ri:
                invalid["complex_too_large_or_empty"] += 1
                continue
            if not (all(g in genes for g in li) and all(g in genes for g in ri)):
                invalid["missing_panel_subunit"] += 1
                continue
            all_pairs.add((li, ri))
    if len(all_pairs) < 20:
        raise ValueError(f"Only {len(all_pairs)} pairs available: check resource name and panel.")
    pairs = sorted(all_pairs)
    ligands = sorted({p[0] for p in pairs})
    receptors = sorted({p[1] for p in pairs})
    needed_genes = sorted({g for li, ri in pairs for g in li + ri})
    l_index = {complex_: i for i, complex_ in enumerate(ligands)}
    r_index = {complex_: i for i, complex_ in enumerate(receptors)}
    mat = csr_matrix((
        np.ones(len(pairs), dtype=np.int16),
        ([l_index[a] for a, b in pairs], [r_index[b] for a, b in pairs])
    ), shape=(len(ligands), len(receptors)), dtype=np.int16)
    mat.sum_duplicates()
    return {
        "pairs": pairs, "ligands": ligands, "receptors": receptors,
        "needed_genes": needed_genes, "adjacency": mat,
        "rows_in_resource": resource_rows, "filter_exclusions": dict(invalid),
    }


def as_fov(x):
    try:
        return int(float(x))
    except (ValueError, TypeError):
        return None


def load_metadata(path, specimens):
    lookup = {}
    allowed = set()
    for geo, cond, lo, hi in specimens:
        allowed.update(range(lo, hi + 1))
    by_fov = defaultdict(dict)
    with path.open(newline="", encoding="utf-8-sig") as fh:
        r = csv.DictReader(fh)
        req = {"fov", "cell", "CenterX_global_px", "CenterY_global_px", "Area"}
        if not req.issubset(set(r.fieldnames or [])):
            raise ValueError(f"Missing metadata fields: {req - set(r.fieldnames or [])}")
        for row in r:
            fov = as_fov(row.get("fov"))
            if fov not in allowed:
                continue
            try:
                x, y, area = (float(row[c]) for c in
                              ("CenterX_global_px", "CenterY_global_px", "Area"))
            except (ValueError, TypeError):
                continue
            if not all(math.isfinite(z) for z in (x, y, area)) or area <= 0:
                continue
            key = str(row["cell"])
            if key in by_fov[fov]:
                raise ValueError(f"Duplicate metadata cell {fov}, {key}")
            by_fov[fov][key] = (x, y, area)
    for geo, _, lo, hi in specimens:
        for f in range(lo, hi + 1):
            lookup[f] = geo
    return by_fov, lookup


def load_counts(path, by_fov, required_genes, detection_count, min_total, chunk_size):
    # Stream chunks of all columns to calculate total RNA count independently
    # of catalogue selection. Keep only binary supported genes per QC cell.
    keep = {f: [] for f in by_fov}
    matched, passing = Counter(), Counter()
    seen = set()
    with pd.read_csv(path, compression="gzip", chunksize=chunk_size,
                     low_memory=False) as reader:
        for ci, chunk in enumerate(reader, 1):
            # Matches columns in the GEO CosMx files, verified in prior audit.
            fovs = pd.to_numeric(chunk["fov"], errors="coerce")
            valid_mask = fovs.isin(by_fov.keys()).to_numpy()
            if not valid_mask.any():
                continue
            subset = chunk.loc[valid_mask].copy()
            if len(subset) == 0:
                continue
            totals = subset.iloc[:, 3:].sum(axis=1, numeric_only=True).to_numpy(dtype=np.float64)
            expr = subset[required_genes].to_numpy(dtype=np.float32, copy=False)
            binary = expr >= detection_count
            s_fovs = fovs.loc[valid_mask].to_numpy(dtype=np.int32)
            s_cells = subset["cell"].astype(str).to_numpy()
            for row_idx, (fov, cell) in enumerate(zip(s_fovs, s_cells)):
                if cell not in by_fov.get(int(fov), {}):
                    continue
                matched[int(fov)] += 1
                if totals[row_idx] < min_total:
                    continue
                signature = (int(fov), cell)
                if signature in seen:
                    raise ValueError(f"Duplicate expression cell {signature}")
                seen.add(signature)
                keep[int(fov)].append((cell, binary[row_idx].copy()))
                passing[int(fov)] += 1
            if ci % 20 == 0:
                print(f"  Read {ci * chunk_size:,} expression rows...", flush=True)
    return keep, matched, passing


def make_complex_matrix(binary, genes, complexes):
    ix = {g: i for i, g in enumerate(genes)}
    cols = []
    for item in complexes:
        if len(item) == 1:
            cols.append(binary[:, ix[item[0]]])
        else:
            cols.append(np.all(binary[:, [ix[g] for g in item]], axis=1))
    # CSR of binary complex-level expression.
    return csr_matrix(np.column_stack(cols).astype(np.int16))


def support_on_pairs(u, v, L, R, A, metric, batch=4000):
    if len(u) == 0:
        return np.empty(0, dtype=np.int32)
    if metric == "ligand":
        # Receptor-active ligand set at every potential receiver.
        Q = R @ A.T
        Q.data[:] = 1
        Q.eliminate_zeros()
        Q = Q.tocsr()
        first, second = L, Q
    else:
        # Pair count: each ligand–receptor complex combination counts once.
        first, second = (L @ A).tocsr(), R
    result = np.empty(len(u), dtype=np.int32)
    for start in range(0, len(u), batch):
        end = min(start + batch, len(u))
        fwd = np.asarray(first[u[start:end]].multiply(second[v[start:end]]).sum(axis=1)).ravel()
        rev = np.asarray(first[v[start:end]].multiply(second[u[start:end]]).sum(axis=1)).ravel()
        result[start:end] = np.maximum(fwd, rev).astype(np.int32)
    return result


def graph_h0(n, edges, support, thresholds):
    # Exact filtration: all vertices enter at time 0; add edges at 1/k
    # where k = supported ligand or LR pair count. Essential components
    # persist to infinity. Arrays contain only positive-support edges.
    order = np.argsort(-support, kind="stable")
    dsu = DSU(n)
    merges = Counter()
    snapshots = []
    idx = 0
    for k in thresholds:
        while idx < len(order) and support[order[idx]] >= k:
            j = int(order[idx])
            if dsu.union(edges[j, 0], edges[j, 1]):
                merges[int(support[j])] += 1
            idx += 1
        snapshots.append({"min_support":int(k), "edges":int(idx),
                          "components":int(dsu.components),
                          "largest_fraction": float(dsu.largest / max(n, 1))})
    return snapshots, merges, dsu


def geometric_graph(n, pairs):
    d = DSU(n)
    for a, b in pairs:
        d.union(a, b)
    return d


def fov_run(coords, signatures, areas, resource, radius_mult, metric,
            rng, permutations, max_neighbors):
    n = len(coords)
    diam = 2 * np.sqrt(areas[areas > 0] / np.pi)
    if not len(diam):
        raise ValueError("No valid cell areas")
    med_diameter = float(np.median(diam))
    radius = med_diameter * radius_mult
    pairs = cKDTree(coords).query_pairs(radius, output_type="ndarray")
    if len(pairs) > max_neighbors:
        raise RuntimeError(f"Too many candidate edges ({len(pairs):,}), limit {max_neighbors:,}")
    if len(pairs):
        pairs = np.sort(pairs.astype(np.int32), axis=1)
        pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]
    else:
        pairs = np.empty((0, 2), dtype=np.int32)
    geo = geometric_graph(n, pairs)
    L = make_complex_matrix(signatures, resource["needed_genes"], resource["ligands"])
    R = make_complex_matrix(signatures, resource["needed_genes"], resource["receptors"])
    u, v = pairs[:, 0], pairs[:, 1]
    strength = support_on_pairs(u, v, L, R, resource["adjacency"], metric)
    valid = strength > 0
    supported_pairs = pairs[valid]
    weights = strength[valid]
    levels = sorted(set(map(int, weights)), reverse=True)
    if 1 not in levels:
        levels.append(1)
    # Every distinct observed strength is included, so no merges are missed.
    curve, mergers, result = graph_h0(n, supported_pairs, weights, levels)
    summary = {
        "cells":n,"spatial_pairs":int(len(pairs)),
        "spatial_component_count":int(geo.components),
        "spatial_largest_fraction":float(geo.largest/max(n,1)),
        "supported_pairs":int(valid.sum()),
        "supported_fraction_of_spatial_pairs":float(valid.mean()) if len(valid) else 0.,
        "positive_support_max":int(weights.max()) if len(weights) else 0,
        "positive_support_median":float(np.median(weights)) if len(weights) else 0.,
        "components_at_k1":int(result.components),
        "largest_fraction_at_k1":float(result.largest/max(n,1)),
        "finite_h0_intervals":int(sum(mergers.values())),
        "essential_h0_intervals":int(result.components),
        "median_cell_diameter_px":med_diameter,
        "spatial_radius_px":float(radius),
    }
    nulls=[]
    for j in range(permutations):
        perm = rng.permutation(n)
        shuf_strength = support_on_pairs(perm[u], perm[v], L, R, resource["adjacency"], metric)
        # Shuffle complete cellular expression profiles relative to measured
        # positions. Preserves cell-cell distances, gene coexpression and
        # per-cell detection, but NOT spatial microenvironment composition.
        good = shuf_strength > 0
        dnull = DSU(n)
        for a,b in pairs[good]:
            dnull.union(a,b)
        nulls.append({"permutation":j+1,
                      "supported_fraction_of_spatial_pairs":float(good.mean()) if len(good) else 0.,
                      "components_at_k1":int(dnull.components),
                      "largest_fraction_at_k1":float(dnull.largest/max(n,1))})
    return summary, curve, mergers, nulls


def write_csv(path, rows, cols):
    with path.open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,
       default=Path.home()/"Real_Data/spatial_ibd_sheaf/pilot_6k")
    p.add_argument("--out",type=Path,
       default=Path.home()/"Real_Data/spatial_ibd_sheaf/outputs/h0_all_fovs")
    p.add_argument("--catalog",type=Path,
       default=Path.home()/"Real_Data/spatial_ibd_sheaf/metadata/liana_omni_resource.csv")
    p.add_argument("--resource", choices=["consensus","cellchatdb","cellphonedb"],
                   default="consensus")
    p.add_argument("--metric", choices=["ligand","pair"],default="ligand")
    p.add_argument("--radius-diameters",type=float,default=2.5)
    p.add_argument("--detection-count",type=float,default=1.)
    p.add_argument("--min-total-count",type=float,default=50.)
    p.add_argument("--permutations",type=int,default=2)
    p.add_argument("--seed",type=int,default=20261031)
    p.add_argument("--chunk-size",type=int,default=384)
    p.add_argument("--max-cells-per-fov",type=int,default=6000)
    p.add_argument("--max-neighbor-pairs",type=int,default=1500000)
    args=p.parse_args()
    if args.radius_diameters<=0 or args.detection_count<=0 or args.permutations<0:
        p.error("radius and detection threshold must be >0, permutations >=0")
    root,out,catalog=args.root.expanduser(),args.out.expanduser(),args.catalog.expanduser()
    out.mkdir(parents=True,exist_ok=True)
    hash_=get_catalog(catalog)
    print("LIANA resource SHA256:",hash_,flush=True)
    all_summaries,all_curves,all_merges,all_nulls,cat_stats=[],[],[],[],[]
    rng=np.random.default_rng(args.seed)

    # Avoid silently merging clinical specimens: each FOV maps to one GEO ID.
    for prefix, samples in SPECIMENS.items():
        expr_path=root/"counts"/(prefix+"_exprMat_file.csv.gz")
        meta_path=root/"metadata"/(prefix+"_metadata_file.csv")
        if not expr_path.exists() or not meta_path.exists():
            raise FileNotFoundError(f"Missing {expr_path} or {meta_path}")
        with pd.read_csv(expr_path, compression="gzip", nrows=0) as hdr:
            gene_columns=list(hdr.columns[3:])
        genes=set(map(str.upper, gene_columns))
        res=catalog_complexes(catalog,genes,args.resource)
        # Actual observed panel gene names are used, preventing case mismatch.
        if len(res["needed_genes"]) < 10:
            raise ValueError("Too few genes from ligand–receptor catalogue")
        if len(set(gene_columns))!=len(gene_columns):
            raise ValueError("Duplicated expression gene columns")
        cat_stats.append({
            "prefix":prefix,"catalogue_entries":res["rows_in_resource"],
            "covered_complex_pairs":len(res["pairs"]),
            "ligand_complexes":len(res["ligands"]),
            "receptor_complexes":len(res["receptors"]),
            "covered_unique_genes":len(res["needed_genes"]),
            "excluded":res["filter_exclusions"]})
        print(f"\n{prefix}: {len(res['pairs']):,} {args.resource} pairs covered; "
              f"{len(res['needed_genes']):,} participating genes.",flush=True)
        meta, fov_to_geo=load_metadata(meta_path,samples)
        print(f"  Loaded coordinates: {sum(len(m) for m in meta.values()):,} cells, "
              f"{len(meta)} FOVs",flush=True)
        records, matched, passing=load_counts(
            expr_path,meta,res["needed_genes"],args.detection_count,
            args.min_total_count,args.chunk_size)
        count_fovs=0
        for fov in sorted(meta):
            geo=fov_to_geo.get(fov)
            if geo is None: continue
            condition=next(s[1] for s in samples if s[0]==geo)
            measured=records.get(fov,[])
            n=len(measured)
            if n<2 or n>args.max_cells_per_fov:
                print(f"  {geo} FOV {fov}: skipped (QC cells={n})",flush=True)
                continue
            coords=np.asarray([meta[fov][cell][:2] for cell,_ in measured],dtype=float)
            areas=np.asarray([meta[fov][cell][2] for cell,_ in measured],dtype=float)
            signatures=np.stack([sig for _,sig in measured])
            stats,curve,merges,nulls=fov_run(
                coords,signatures,areas,res,args.radius_diameters,
                args.metric,rng,args.permutations,args.max_neighbor_pairs)
            extra={"geo":geo,"condition":condition,"prefix":prefix,"fov":int(fov),
                   "metadata_cells":len(meta[fov]),"matched_expression_cells":int(matched[fov]),
                   "passing_qc_expression_cells":int(passing[fov])}
            all_summaries.append({**extra,**stats})
            all_curves.extend({**extra,"min_support":k["min_support"],
                               "filtration_weight":1.0/k["min_support"],
                               "edges":k["edges"],"components":k["components"],
                               "largest_fraction":k["largest_fraction"]}
                              for k in curve)
            all_merges.extend({**extra,"support":int(k),
                               "death":1.0/int(k),
                               "n_merges":int(v)}
                              for k,v in sorted(merges.items(),reverse=True))
            all_nulls.extend({**extra,**null} for null in nulls)
            count_fovs+=1
            print(f"  {geo}, FOV {fov}: n={n}, nearby={stats['spatial_pairs']:,}, "
                  f"LR={stats['supported_pairs']:,}, "
                  f"giant={stats['largest_fraction_at_k1']:.3f}, "
                  f"components={stats['components_at_k1']}",flush=True)
        print(f"Finished {prefix}: {count_fovs} FOVs",flush=True)

    if not all_summaries:
        raise RuntimeError("No eligible FOVs; inspect QC and file joins.")
    write_csv(out/"fov_summary.csv",all_summaries,list(all_summaries[0]))
    write_csv(out/"component_filtration.csv",all_curves,list(all_curves[0]))
    write_csv(out/"finite_h0_death_histogram.csv",all_merges,list(all_merges[0]))
    if all_nulls:
        write_csv(out/"expression_shuffle_nulls.csv",all_nulls,list(all_nulls[0]))
    grouped=defaultdict(list)
    for row in all_summaries:
        grouped[row["geo"]].append(row)
    aggregates=[]
    for geo,rows in grouped.items():
        # Descriptive FOV aggregation only; NOT inferential replicate counts.
        total_nbr=sum(r["spatial_pairs"] for r in rows)
        aggregates.append({
            "geo":geo,"condition":rows[0]["condition"],
            "n_fovs":len(rows),"n_cells":sum(r["cells"] for r in rows),
            "spatial_pairs":total_nbr,
            "supported_pairs":sum(r["supported_pairs"] for r in rows),
            "supported_fraction_pooled":sum(r["supported_pairs"] for r in rows)/max(1,total_nbr),
            "median_fov_largest_fraction":float(median(
                r["largest_fraction_at_k1"] for r in rows)),
            "median_fov_geometric_largest_fraction":float(median(
                r["spatial_largest_fraction"] for r in rows)),
        })
    write_csv(out/"specimen_summary.csv",aggregates,list(aggregates[0]))
    (out/"run_config.json").write_text(json.dumps({
        "parameters":{k:str(v) if isinstance(v,Path) else v
                       for k,v in vars(args).items()},
        "catalogue_url":CATALOG_URL,"catalogue_sha256":hash_,
        "catalogue_filters":cat_stats,
        "interpretation":"Exploratory candidate communication connectivity, not actual transport, sheaf compatibility or response."
           "Complex subunits must all be present. Edges eligible by FOV-centroid distance and gene detection."
           "Undirected H0 OR of directional evidence, weight 1/(distinct ligand complexes or LR complexes)."
           "Each sampled FOV is nested inside one of only four GEO specimens; FOVs are not independent donors."
           "Null shuffles complete binary expression profiles among fixed within-FOV spatial positions."
           "The catalogue can contain interactions beyond canonical secreted signalling and does not establish activity."
    },indent=2),encoding="utf-8")
    print("\nFinal specimen-level descriptive summary:",flush=True)
    for s in aggregates:
        print(f" {s['geo']} {s['condition']}: {s['n_fovs']} FOVs, "
              f"{s['n_cells']:,} cells, supported spatial pairs "
              f"{s['supported_fraction_pooled']:.2%}, "
              f"median FOV giant {s['median_fov_largest_fraction']:.3f}",flush=True)
    print(f"Files saved under {out}",flush=True)

if __name__=="__main__":
    main()
