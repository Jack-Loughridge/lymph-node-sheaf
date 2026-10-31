#!/usr/bin/env python3
"""Audits Capeling et al. Supplementary Data 7 response genes against CosMx 6K.

Paper: https://doi.org/10.1038/s41467-025-68247-6
Supplementary Data 7: cell-type-specific DEGs after factor stimulation.

Uses Python standard-library XLSX streaming (zipfile+ElementTree).
Pandas/NumPy/openpyxl are NOT required. Never exports underlying cell data.

Run:
  python3 scripts/audit_organoid_response_genes.py
  python3 scripts/audit_organoid_response_genes.py --supplement /path/to/SuppData7.xlsx

Caveat: 'responsive factors' are perturbations in inflamed epithelial organoids;
these are NOT automatically active in the CosMx tissue. This audit measures
TARGET GENE PANEL COVERAGE, not dose-response calibration or actual signalling.
"""
import argparse
import csv
from collections import defaultdict
import json
import math
from pathlib import Path
import posixpath
import re
import shutil
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile, BadZipFile

DOI = "10.1038/s41467-025-68247-6"
# Nature supplementary numbering: information=1, review=2,
# additional-file legend=3, Data 1=4, hence Data 7=10.
SUPP_URL = ("https://static-content.springer.com/esm/"
            "art%3A10.1038%2Fs41467-025-68247-6/MediaObjects/"
            "41467_2025_68247_MOESM10_ESM.xlsx")
FACTOR_FILE = "79_factor_panel_and_liana_coverage.csv"
ALIASES = {"TNFA":"TNF", "ISlR":"ISLR", "IL1ALPHA":"IL1A",
           "TWEAK":"TNFSF12", "TRAIL":"TNFSF10", "TNFALPHA":"TNFA",
           "IFNGAMMA":"IFNG", "CT1":"CTF1", "NT4":"NTF4"}
NS = {"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "p":"http://schemas.openxmlformats.org/package/2006/relationships"}
BASE = "{"+NS["m"]+"}"
GENE_COLS = ["gene","genes","gene_symbol","genesymbol","feature","symbol"]
PERT_COLS = ["pert","perturbation","treatment","factor","ligand","condition"]
P_COLS = ["padj","p_val_adj","pvaladj","adjustedpvalue","fdr","p_adj",
          "pvalueadjusted","p.adj"]
FC_COLS = ["log2foldchange","avg_log2fc","avglog2fc","log2fc","logfc",
           "avg_logfc","log2_fold_change"]
CELL_COLS = ["celltype","cell_type","cell.type","cell_type_annotation",
             "cluster","celltypeannotation","cell","cell_type_name"]
OUTPUT_FIELDS = [
    "factor","ligand_gene_in_panel","complete_receptor_in_panel",
    "significant_fdr_genes","significant_fdr_genes_in_panel",
    "significant_fdr_abs_log2fc_gt_1_genes",
    "response_genes_in_panel","response_coverage_fraction",
    "cell_types_with_significant_response","cell_types_with_5_covered_response_genes",
    "example_covered_response_genes"
]

def normalize(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())

def factor_name(value):
    text=str(value or "").strip()
    text=re.sub(r"[\s_\-]+","",text).upper().replace("Α","A").replace("Γ","G")
    canonical={normalize(k):v for k,v in ALIASES.items()}
    return canonical.get(normalize(text), text)

def download_supplement(path, url):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.is_file() and path.stat().st_size > 10_000_000:
        print(f"Using existing Supplementary Data 7: {path}",flush=True)
        return
    temp=path.with_suffix(".part")
    for trial in range(3):
        try:
            req=Request(url,headers={
                "User-Agent":"Mozilla/5.0 (compatible; research-audit/1.0)",
                "Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*"})
            with urlopen(req,timeout=90) as response, temp.open("wb") as out:
                while True:
                    block=response.read(1024*1024)
                    if not block:break
                    out.write(block)
                    if out.tell()>160_000_000:
                        raise RuntimeError("Unexpectedly large supplement (>160MB)")
            with ZipFile(temp) as z:
                if "xl/workbook.xml" not in z.namelist():
                    raise ValueError("Download was not the expected XLSX workbook")
            temp.replace(path)
            print(f"Downloaded {path} ({path.stat().st_size/1024**2:.1f} MiB)")
            return
        except (HTTPError,URLError,ValueError,BadZipFile,RuntimeError) as exc:
            temp.unlink(missing_ok=True)
            print(f"Supplement download attempt {trial+1}: {exc}",flush=True)
            time.sleep(trial+1)
    raise RuntimeError(
        "Cannot automatically download Supplementary Data 7. Download it "
        "from the Nature article's Supplementary Data 7 link, then rerun with "
        "--supplement /path/to/file.xlsx. No data have been deleted.")

def string_table(z):
    name="xl/sharedStrings.xml"
    if name not in z.namelist():return []
    result=[]
    with z.open(name) as fh:
        for _,element in ET.iterparse(fh,events=("end",)):
            if element.tag==BASE+"si":
                value="".join((t.text or "") for t in element.iter(BASE+"t"))
                result.append(value)
                element.clear()
    return result

def excel_col(ref):
    chars=re.match(r"([A-Z]+)",ref or "")
    if chars is None:return None
    n=0
    for letter in chars.group(1):
        n=n*26+ord(letter)-64
    return n-1

def sheets(z):
    workbook=ET.fromstring(z.read("xl/workbook.xml"))
    rels=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    relmap={r.attrib["Id"]:r.attrib.get("Target","") for r in rels}
    for entry in workbook.findall("m:sheets/m:sheet", NS):
        name=entry.attrib["name"]
        rid=entry.attrib.get("{"+NS["r"]+"}id")
        target=relmap.get(rid,"")
        if target.startswith("/"): target=target.lstrip("/")
        elif not target.startswith("xl/"):
            target=posixpath.normpath(posixpath.join("xl",target))
        if target in z.namelist():
            yield name,target

def cells_by_row(z,path,strings):
    with z.open(path) as fh:
        context=ET.iterparse(fh,events=("start","end"))
        root=None
        for event,element in context:
            if root is None and event=="start":root=element
            if event!="end" or element.tag!=BASE+"row":continue
            values={}
            for cell in element.findall("m:c",NS):
                ix=excel_col(cell.attrib.get("r",""))
                if ix is None:continue
                typ=cell.attrib.get("t")
                v=cell.find("m:v",NS)
                inline=cell.find("m:is",NS)
                if typ=="s" and v is not None:
                    try:raw=strings[int(v.text)]
                    except (ValueError,IndexError,TypeError):raw=""
                elif typ=="inlineStr" and inline is not None:
                    raw="".join(t.text or "" for t in inline.iter(BASE+"t"))
                elif v is not None:
                    raw=v.text or ""
                else:
                    raw=""
                values[ix]=raw
            yield values
            element.clear()
            # Drops accumulated cleared row nodes to bound memory.
            if root is not None:root.clear()

def find_header(row):
    columns={normalize(text):idx for idx,text in row.items()}
    def pick(aliases):
        return next((columns[normalize(a)] for a in aliases
                     if normalize(a) in columns),None)
    return {"gene":pick(GENE_COLS),"pert":pick(PERT_COLS),
            "padj":pick(P_COLS),"logfc":pick(FC_COLS),
            "celltype":pick(CELL_COLS)}

def fnumber(value):
    try:return float(str(value))
    except (ValueError,TypeError):return float("nan")

def read_panel(root):
    files=list((root/"counts").glob("*_exprMat_file.csv.gz"))
    if not files:raise FileNotFoundError(
        f"Cannot find CosMx expression headers in {root/'counts'}")
    import gzip
    panels=[]
    for path in files:
        with gzip.open(path,"rt",newline="",encoding="utf-8-sig") as f:
            header=next(csv.reader(f))
        panels.append({gene.upper() for gene in header[3:]})
    return set.intersection(*panels),[str(f) for f in files]

def read_factor_metadata(path):
    if not path.exists():
        raise FileNotFoundError(f"Need previous 79-factor audit table: {path}")
    with path.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def writecsv(path,rows,headers):
    with path.open("w",newline="",encoding="utf-8") as fh:
        wr=csv.DictWriter(fh,fieldnames=headers,extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)

def audit(args):
    out=args.out.expanduser()
    out.mkdir(parents=True,exist_ok=True)
    panel,sourcefiles=read_panel(args.root.expanduser())
    prior=read_factor_metadata(args.previous.expanduser())
    factors={factor_name(r["factor_reagent_name"]):r for r in prior}
    # Protein reagents include ambiguous IFNA, so preserve as a treatment label.
    source=args.supplement.expanduser()
    if not args.no_download:
        download_supplement(source,args.url)
    if not source.exists():
        raise FileNotFoundError(source)
    measured=defaultdict(lambda:{"fdr":set(),"stringent":set(),"cells":defaultdict(set)})
    sheet_summary=[]
    examples=defaultdict(dict)
    with ZipFile(source) as z:
        strings=string_table(z)
        workbook_sheets=list(sheets(z))
        print(f"Supplement contains {len(workbook_sheets)} worksheets; "
              f"{len(strings)} shared strings.",flush=True)
        for name,xmlpath in workbook_sheets:
            headers=None
            processed=0
            sigcount=0
            preview=[]
            recognized=False
            for index,row in enumerate(cells_by_row(z,xmlpath,strings)):
                if index<4:preview.append(" | ".join(str(x) for x in list(row.values())[:9])[:350])
                if headers is None:
                    h=find_header(row)
                    if h["gene"] is not None and h["padj"] is not None:
                        headers=h
                        recognized=True
                        continue
                    if index>=25:break
                    continue
                gene=str(row.get(headers["gene"],"")).strip().upper()
                if not re.fullmatch(r"[A-Z0-9.-]{2,35}",gene):continue
                pertraw=row.get(headers["pert"],"") if headers["pert"] is not None else name
                factor=factor_name(pertraw)
                if factor not in factors:
                    # If the worksheet is named for one experimental factor, use it.
                    fac2=factor_name(name)
                    if fac2 in factors:factor=fac2
                    else:continue
                padj=fnumber(row.get(headers["padj"],""))
                fc=fnumber(row.get(headers["logfc"],"")) if headers["logfc"] is not None else float("nan")
                celltype=str(row.get(headers["celltype"],"")).strip() if headers["celltype"] is not None else name
                processed+=1
                if math.isfinite(padj) and padj<.05:
                    measured[factor]["fdr"].add(gene)
                    if math.isfinite(fc) and abs(fc)>1:
                        measured[factor]["stringent"].add(gene)
                        measured[factor]["cells"][celltype].add(gene)
                        if gene in panel:
                            examples[factor][gene]=max(abs(fc),examples[factor].get(gene,0))
                        sigcount+=1
            sheet_summary.append({
                "sheet":name,"headers_found":recognized,
                "pert_column":headers["pert"] if headers else None,
                "gene_column":headers["gene"] if headers else None,
                "padj_column":headers["padj"] if headers else None,
                "logfc_column":headers["logfc"] if headers else None,
                "accepted_factor_gene_rows":processed,
                "significant_factor_gene_rows":sigcount,
                "preview":" || ".join(preview)[:500]})
            print(f"{name}: parsed={processed:,}; significant={sigcount:,}; "
                  f"headers={'yes' if recognized else 'NO'}",flush=True)
    rows=[]
    bycell=[]
    for factor,info in sorted(factors.items()):
        record=measured[factor]
        strong=record["stringent"]
        sig=record["fdr"]
        included=strong & panel
        celltypecount=0
        for ct, genes in sorted(record["cells"].items()):
            overlap=genes&panel
            if len(overlap)>=5:celltypecount+=1
            bycell.append({
                "factor":factor,"receiver_cell_type":ct,
                "n_response_genes":len(genes),
                "n_response_genes_measured":len(overlap),
                "coverage":len(overlap)/len(genes) if genes else 0,
                "top_covered_genes":";".join(sorted(overlap)[:30])
            })
        rows.append({
            "factor":factor,
            "ligand_gene_in_panel":info.get("factor_gene_in_cosmx_panel",""),
            "complete_receptor_in_panel":info.get("full_panel_lr_eligible",""),
            "significant_fdr_genes":len(sig),
            "significant_fdr_genes_in_panel":len(sig&panel),
            "significant_fdr_abs_log2fc_gt_1_genes":len(strong),
            "response_genes_in_panel":len(included),
            "response_coverage_fraction":round(len(included)/len(strong),4) if strong else 0,
            "cell_types_with_significant_response":len(record["cells"]),
            "cell_types_with_5_covered_response_genes":celltypecount,
            "example_covered_response_genes":";".join(
                gene for gene,weight in sorted(examples[factor].items(),
                   key=lambda x:-x[1])[:25])
        })
    writecsv(out/"response_overlap_by_factor.csv",rows,OUTPUT_FIELDS)
    writecsv(out/"response_overlap_by_celltype.csv",bycell,
       ["factor","receiver_cell_type","n_response_genes",
        "n_response_genes_measured","coverage","top_covered_genes"])
    writecsv(out/"supplement_sheet_audit.csv",sheet_summary,
       list(sheet_summary[0]) if sheet_summary else ["sheet"])
    shortlist=[x for x in rows
        if int(x["complete_receptor_in_panel"] or 0)>0
        and x["response_genes_in_panel"]>=5
        and x["cell_types_with_5_covered_response_genes"]>=1]
    writecsv(out/"provisional_compatibility_candidates.csv",shortlist,OUTPUT_FIELDS)
    summary={
        "measured_genes_in_CosMx_6K":len(panel),
        "factor_treatments_audited":len(factors),
        "factors_with_any_fdr_gene":sum(x["significant_fdr_genes"]>0 for x in rows),
        "factors_with_strong_response":sum(
           x["significant_fdr_abs_log2fc_gt_1_genes"]>=2 for x in rows),
        "factors_with_5_measured_strong_response_genes":sum(
           x["response_genes_in_panel"]>=5 for x in rows),
        "factors_with_5_measured_response_genes_and_complete_receptor":len(shortlist),
        "sheets_usable":sum(x["accepted_factor_gene_rows"]>0 for x in sheet_summary),
        "sheets_inspected":len(sheet_summary),
        "source":{"paper_doi":DOI,"supplement_url":args.url,
                  "local_supplement":str(source),
                  "CosMx_headers":sourcefiles},
        "thresholds":{"padj_lt":0.05,"abs_log2fc_gt":1,
                      "minimum_measured_targets_for_candidate":5},
        "caveats":[
           "Panel overlap is not gene detection or receiver predictive validity.",
           "The 79 perturbations are extracellular additions in three epithelial organoid donors; not direct signalling measurements.",
           "Effects are context-dependent (cytomix) and high-dose for many factors.",
           "Shortlist criterion is a feasibility screen, NOT independently validated compatibility.",
           "Some factor name aliases or supplementary formats may need manual curation.",
           "Spatial allocation by COMMOT, receptor availability and validation still outstanding."
        ]}
    (out/"response_overlap_summary.json").write_text(json.dumps(summary,indent=2))
    print("\nFINAL SUMMARY\n"+json.dumps(summary,indent=2),flush=True)
    if summary["sheets_usable"]==0:
        raise RuntimeError("No suitable tables parsed. Inspect supplement_sheet_audit.csv "
                           "and adjust column aliases before interpreting zero counts.")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,
       default=Path.home()/"Real_Data/spatial_ibd_sheaf/pilot_6k")
    p.add_argument("--previous",type=Path,
       default=Path.home()/"Real_Data/spatial_ibd_sheaf/outputs/organoid_79_audit/"+FACTOR_FILE)
    p.add_argument("--out",type=Path,
       default=Path.home()/"Real_Data/spatial_ibd_sheaf/outputs/organoid_response_audit")
    p.add_argument("--supplement",type=Path,
       default=Path.home()/"Real_Data/spatial_ibd_sheaf/metadata/Capeling_Supplementary_Data_7.xlsx")
    p.add_argument("--url",default=SUPP_URL)
    p.add_argument("--no-download",action="store_true")
    args=p.parse_args()
    audit(args)

if __name__=="__main__":
    main()
