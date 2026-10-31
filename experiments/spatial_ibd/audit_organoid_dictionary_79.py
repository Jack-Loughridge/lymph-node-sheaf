#!/usr/bin/env python3
"""Audit Capeling 79-factor intestinal organoid dictionary against CosMx 6K.

Sources:
  Capeling et al., Nature Communications 2026, DOI 10.1038/s41467-025-68247-6
  Factor names extracted from Genentech/secretome_dictionary,
  Cytomix_Secretome_Processing.Rmd (79 distinct ligand treatments).
  LIANA consensus ligand-receptor resource, same cached file as prior H0 run.

Separately reports (1) ligand genes present on panel, (2) complete receptor
complexes captured by panel, (3) observed ligand and receptor expression in
each of four specimens, and (4) fraction of LIANA catalogue interaction pairs
covered by the experimental factor subset.

Does NOT infer protein secretion, physical ligand concentration, dose-response,
epithelial receptor function, or whether downstream perturbation DEGs are
covered by the CosMx panel. These require separate analysis.

Output contains aggregate data only. No cell-level data are exported.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

RAW_FACTORS = (
    "ANGPTL1 ANGPTL4 AREG BDNF BMP2 BMP3 BMP4 BMP5 BMP6 BMP8B BTC "
    "CTF1 CTSG DKK1 EDIL3 EDN1 EGF EREG F2 FGF1 FGF18 FGF2 FGF7 "
    "FGF9 GDF10 GDF11 GDF15 GDF5 GDF6 GH1 GRN GUCA2A HBEGF HGF "
    "IFNA IFNE IFNG IGF1 IGF2 IHH IL13 IL17B IL1A IL22 IL4 INHBA "
    "INS ISlR KISS1 KITLG LGALS9 LIF MDK MIF NAMPT NPPC NRG1 "
    "NTF3 NTF4 OSM POSTN PPY PRL PROS1 PTN SEMA3C SHH SLIT3 SPP1 "
    "SST TGFA TGFB1 TGFB2 TGFB3 TNFA TNFSF10 TNFSF11 TNFSF12 VIP"
).split()

assert len(RAW_FACTORS)==79 and len(set(RAW_FACTORS))==79
# Authors' reagent names, not always HGNC symbols.
ALIASES={"TNFA":"TNF","ISlR":"ISLR"}
AMBIGUOUS={"IFNA"} # IFNA denotes the interferon-alpha family, not one HGNC gene

SPECIMENS={
    "Cos1400R355RRNA6K4u1716":[
        ("GSM9345416","healthy",32,65),
        ("GSM9345417","uc_noninflamed",1,31)],
    "Cos1380R387RRNA6K4u2716":[
        ("GSM9345418","uc_inflamed",18,45),
        ("GSM9345419","uc_inflamed",1,17)],
}
def write_csv(path, rows):
    if rows:
        with path.open("w", newline="", encoding="utf-8") as fh:
            w=csv.DictWriter(fh, fieldnames=list(rows[0]),extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

def receptor_map(resource_path, genes):
    """Conservative match of ligand single-gene names to LIANA receptors."""
    out=defaultdict(set)
    all_covered_pairs=set()
    with resource_path.open(newline="",encoding="utf-8-sig") as fh:
        rd=csv.DictReader(fh)
        for row in rd:
            if str(row.get("resource","")).lower().strip()!="consensus": continue
            ligand=str(row.get("source_genesymbol","")).strip().upper()
            receptor=str(row.get("target_genesymbol","")).strip().upper()
            if not ligand or not receptor: continue
            subunits=tuple(receptor.split("_"))
            if not all(re.fullmatch(r"[A-Z0-9-]+", x) for x in subunits):
                continue
            # Preserve receptor complexes; no partial receptor matching.
            out[ligand].add(subunits)
            if (ligand in genes and all(g in genes for g in subunits)):
                all_covered_pairs.add((ligand,subunits))
    return out,all_covered_pairs

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/pilot_6k")
    p.add_argument("--catalog",type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/metadata/liana_omni_resource.csv")
    p.add_argument("--out",type=Path,
        default=Path.home()/"Real_Data/spatial_ibd_sheaf/outputs/organoid_79_audit")
    p.add_argument("--min-total-umis",type=int,default=50)
    p.add_argument("--chunksize",type=int,default=256)
    p.add_argument("--panel-only",action="store_true",
                   help="Report measured ligand and receptor genes without scanning all cell counts")
    args=p.parse_args()
    root=args.root.expanduser(); out=args.out.expanduser()
    out.mkdir(parents=True,exist_ok=True)
    matrixpaths=[
        root/"counts"/(prefix+"_exprMat_file.csv.gz")
        for prefix in SPECIMENS
    ]
    gene_panels=[]
    for path in matrixpaths:
        if not path.is_file():raise FileNotFoundError(path)
        with gzip.open(path,"rt",newline="",encoding="utf-8-sig") as f:
            header=next(csv.reader(f))
        if header[:3]!=["fov","cell_ID","cell"]:
            raise RuntimeError("Unexpected gene table columns: "+str(header[:3]))
        gene_panels.append(set(g.upper() for g in header[3:]))
    if gene_panels[0]!=gene_panels[1]:
        raise RuntimeError("The two 6K specimens have different gene panels")
    panel=gene_panels[0]
    if not args.catalog.is_file():
        raise FileNotFoundError("Missing LIANA resource: "+str(args.catalog))
    receptors, all_catalogue_pairs=receptor_map(args.catalog,panel)

    rows=[]
    for name in RAW_FACTORS:
        if name in AMBIGUOUS:
            symbol=None; status="ambiguous_factor_family"
        else:
            symbol=ALIASES.get(name,name).upper()
            status="single_gene_symbol"
        all_rec=sorted(receptors.get(symbol,set())) if symbol else []
        covered_rec=[r for r in all_rec if all(g in panel for g in r)]
        rows.append({
            "factor_reagent_name":name,
            "gene_symbol":symbol or "",
            "status":status,
            "factor_gene_in_cosmx_panel":int(bool(symbol and symbol in panel)),
            "receptor_complexes_listed_in_liana":len(all_rec),
            "complete_receptor_complexes_in_cosmx_panel":len(covered_rec),
            "full_panel_lr_eligible":int(bool(symbol and symbol in panel and covered_rec)),
            "example_panel_receptor_complexes":";".join("_".join(r) for r in covered_rec[:15]),
        })
    if args.panel_only:
        write_csv(out/"79_factor_panel_and_liana_coverage.csv",rows)
        ligand_measured=sum(x["factor_gene_in_cosmx_panel"] for x in rows)
        LR_eligible=sum(x["full_panel_lr_eligible"] for x in rows)
        missing=[x["factor_reagent_name"] for x in rows
                 if x["status"]=="single_gene_symbol"
                 and not x["factor_gene_in_cosmx_panel"]]
        summary={
            "number_of_treatments":79,
            "ambiguous_gene_family":"IFNA (cannot represent by one HGNC symbol)",
            "single_gene_factors":78,
            "ligand_genes_in_CosMx_panel":ligand_measured,
            "factors_with_ligand_and_complete_LIANA_receptor_in_panel":LR_eligible,
            "not_measured_ligand_genes":missing,
            "note":"Panel presence only: does not establish actual expression, protein secretion, or downstream response."
        }
        (out/"panel_only_summary.json").write_text(json.dumps(summary,indent=2))
        print(json.dumps(summary,indent=2))
        print("See",out/"79_factor_panel_and_liana_coverage.csv")
        return
    flookup={x["factor_reagent_name"]:x for x in rows}
    sample_counts={}
    stats={ (geo,name):[0,0] for group in SPECIMENS.values()
            for geo,_,_,_ in group for name in RAW_FACTORS}
    total_qc=Counter()
    sampled={}

    for path, (prefix, spec) in zip(matrixpaths,SPECIMENS.items()):
        fov_to_geo={f:geo for geo,cond,lo,hi in spec for f in range(lo,hi+1)}
        allowed=set(fov_to_geo)
        # Load all genes only in small chunks to compute the same QC threshold
        # as our H0 pilot, and then count factor/receptor expression.
        want={s for x in rows if x["status"]=="single_gene_symbol"
              for s in [x["gene_symbol"]] if s in panel}
        for x in rows:
            if x["gene_symbol"]:
                for receptor in receptors.get(x["gene_symbol"],()):
                    if all(g in panel for g in receptor):
                        want.update(receptor)
        want=sorted(want)
        print(f"Reading {path.name}: {len(want)} factor/receptor genes.",flush=True)
        for chunk in pd.read_csv(path,compression="gzip",chunksize=args.chunksize,
                                 low_memory=False):
            fov=pd.to_numeric(chunk["fov"],errors="coerce")
            selected=fov.isin(allowed).to_numpy()
            if not selected.any(): continue
            sub=chunk.loc[selected]
            samples=fov.loc[selected].astype(int).map(fov_to_geo).to_numpy()
            total=sub.iloc[:,3:].sum(axis=1,numeric_only=True).to_numpy()
            passed=(total>=args.min_total_umis)
            data=sub.loc[:,want].to_numpy(dtype=np.float32,copy=False)>0
            ix={g:i for i,g in enumerate(want)}
            positive={g:data[:,i] for g,i in ix.items()}
            for geo,_,_,_ in spec:
                selector=(samples==geo)&passed
                if not selector.any(): continue
                count=int(selector.sum())
                total_qc[geo]+=count
                for name in RAW_FACTORS:
                    rec=flookup[name]
                    sym=rec["gene_symbol"]
                    ligand_n=int(positive[sym][selector].sum()) if sym in positive else 0
                    receiver=np.zeros(len(data),dtype=bool)
                    for complex_ in receptors.get(sym,()):
                        if not all(g in panel for g in complex_):continue
                        condition=np.ones(len(data),dtype=bool)
                        for g in complex_:
                            condition &= positive[g]
                        receiver |= condition
                    receiver_n=int(receiver[selector].sum())
                    old=stats[(geo,name)]
                    old[0]+=ligand_n;old[1]+=receiver_n
        print(f"  QC cell counts: "+", ".join(f"{geo}:{total_qc[geo]}"
                                           for geo,_,_,_ in spec),flush=True)

    bysample=[]
    for prefix, spec in SPECIMENS.items():
        for geo,condition,_,_ in spec:
            n=total_qc[geo]
            for name in RAW_FACTORS:
                factor_n, receptor_n=stats[(geo,name)]
                bysample.append({
                    "geo":geo,"condition":condition,
                    "factor":name,"n_qc_cells":n,
                    "sender_factor_gene_positive_cells":factor_n,
                    "sender_gene_detection_fraction":factor_n/max(n,1),
                    "any_complete_receptor_positive_cells":receptor_n,
                    "receptor_detection_fraction":receptor_n/max(n,1)
                })
    for row in rows:
        all_counts=[stats[(geo,row["factor_reagent_name"])]
                    for group in SPECIMENS.values() for geo,_,_,_ in group]
        row["factor_gene_detected_in_any_qc_cell"]=int(any(x[0] for x in all_counts))
        row["receptor_detected_in_any_qc_cell"]=int(any(x[1] for x in all_counts))
        row["both_detected_in_cohort"]=int(
            row["full_panel_lr_eligible"] and any(x[0] for x in all_counts)
            and any(x[1] for x in all_counts)
        )
        row["total_ligand_positive_qc_cells"]=sum(x[0] for x in all_counts)

    direct=[r for r in rows if r["status"]=="single_gene_symbol"]
    covered=[r for r in rows if r["full_panel_lr_eligible"]]
    eligible={r["gene_symbol"] for r in covered}
    pairs_in_dictionary=[p for p in all_catalogue_pairs if p[0] in eligible]
    totals={
        "n_factor_treatments":79,
        "n_unambiguous_single_gene_symbols":len(direct),
        "factor_family_ambiguous":[r["factor_reagent_name"] for r in rows
                                   if r["status"]!="single_gene_symbol"],
        "ligand_genes_measured_count":sum(r["factor_gene_in_cosmx_panel"] for r in direct),
        "ligand_genes_measured_denominator":len(direct),
        "factor_with_receptor_complex_covered":len(covered),
        "factor_ligand_detected_in_any_cosmx_qc_cell":sum(
            r["factor_gene_detected_in_any_qc_cell"] for r in direct),
        "factor_with_ligand_and_receptor_detected_in_cohort":sum(
            r["both_detected_in_cohort"] for r in rows),
        "n_distinct_liana_pairs_on_panel_for_dictionary_ligands":len(pairs_in_dictionary),
        "n_liana_single_ligand_pairs_on_panel":len(all_catalogue_pairs),
        "share_liana_single_ligand_pairs_from_dictionary":
            len(pairs_in_dictionary)/max(1,len(all_catalogue_pairs)),
        "qc_cells_by_sample":dict(total_qc),
        "critical_limitations":[
            "Ligand gene on panel does not mean ligand protein secreted.",
            "Cytokine stimulation response is measured only for intestinal epithelial organoid receiver cells.",
            "LIANA coverage is only one curated resource; absence does not prove no receptor.",
            "Receptor-gene detection does not establish receptor protein or active pathway.",
            "Measured expression and receiver responses are not simultaneously observed in the same organoid-spatial cells.",
            "No spatial allocation or COMMOT calculation is performed by this audit.",
            "Perturbation-specific response target coverage requires supplementary DEG / gene-program tables, not just factor genes.",
        ]
    }
    write_csv(out/"79_factor_panel_and_liana_coverage.csv",rows)
    write_csv(out/"79_factor_expression_by_specimen.csv",bysample)
    (out/"summary.json").write_text(json.dumps(totals,indent=2),encoding="utf-8")
    print("\nFACTOR COVERAGE RESULTS")
    print("79 factors, of which 1 (IFNA) is a non-unique family designation.")
    for key in ("ligand_genes_measured_count",
                "factor_with_receptor_complex_covered",
                "factor_ligand_detected_in_any_cosmx_qc_cell",
                "factor_with_ligand_and_receptor_detected_in_cohort",
                "n_distinct_liana_pairs_on_panel_for_dictionary_ligands",
                "n_liana_single_ligand_pairs_on_panel"):
        print(f"{key}: {totals[key]}")
    print("\nMissing gene targets:",", ".join(r["factor_reagent_name"]
        for r in direct if not r["factor_gene_in_cosmx_panel"]) or "none")
    print("\nWrote:",out)
    print("Do not equate this coverage with 79 validated functional communication pathways.")


if __name__=="__main__":
    main()
