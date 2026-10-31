#!/usr/bin/env python3
"""Audit signed Capeling organoid DE responses against existing CosMx counts.

Inputs stay read-only. Export aggregate evidence only, to a fresh output directory.
Both value-identical worksheets are quarantined from conservative candidates.
Panel coverage, observed transcript detection, and signalling activity are distinct.
Requires numpy/pandas for counts; tests/validation additionally use openpyxl.
Paper: https://doi.org/10.1038/s41467-025-68247-6
"""
import argparse
import csv
from collections import defaultdict
import json
import hashlib
from collections import Counter
import math
from pathlib import Path
import posixpath
import re
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
           "TWEAK":"TNFSF12", "TRAIL":"TNFSF10", "TNFALPHA":"TNF",
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

def normalize(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())

CANONICAL_ALIASES = {normalize(k):v for k,v in ALIASES.items()}

def factor_name(value):
    text=str(value or "").strip()
    text=re.sub(r"[\s_\-]+","",text).upper().replace("Α","A").replace("Γ","G")
    return CANONICAL_ALIASES.get(normalize(text), text)

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

def technical_target(gene):
    return bool(re.fullmatch(r"(?:NEGATIVE|SYSTEMCONTROL|NEGPRB|NEGPROBE|FALSECODE)[0-9_-]*", gene.upper()))

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
    if any(p != panels[0] for p in panels[1:]):
        raise ValueError("CosMx panels differ; audit each panel separately")
    controls = {g for g in panels[0] if technical_target(g)}
    return panels[0] - controls,[str(f) for f in sorted(files)]

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

# Conservative signalling receptors, not every binding partner in LIANA.
# See receptor_sources.json emitted with the audit. Alternatives are OR; subunits AND.
CANONICAL_RECEPTORS = {
    'BMP2': [('BMPR1A','BMPR2'), ('BMPR1B','BMPR2'),
             ('BMPR1A','ACVR2A'), ('BMPR1A','ACVR2B'),
             ('BMPR1B','ACVR2A'), ('BMPR1B','ACVR2B')],
    'BMP4': [('BMPR1A','BMPR2'), ('BMPR1B','BMPR2'),
             ('BMPR1A','ACVR2A'), ('BMPR1A','ACVR2B'),
             ('BMPR1B','ACVR2A'), ('BMPR1B','ACVR2B')],
    'BMP6': [('ACVR1','ACVR2A'), ('ACVR1','ACVR2B'), ('ACVR1','BMPR2')],
    'IFNG': [('IFNGR1','IFNGR2')],
    'IL13': [('IL13RA1','IL4R')],
    # Type II receptor is the relevant epithelial choice. Type I is separately noted.
    'IL4': [('IL13RA1','IL4R')],
    'IL22': [('IL22RA1','IL10RB')],
    'INHBA': [('ACVR1B','ACVR2A'), ('ACVR1B','ACVR2B')],
    'TGFB1': [('TGFBR1','TGFBR2')],
    'TGFB3': [('TGFBR1','TGFBR2')],
    'TNF': [('TNFRSF1A',), ('TNFRSF1B',)],
    'TNFSF10': [('TNFRSF10A',), ('TNFRSF10B',)],
    'TNFSF11': [('TNFRSF11A',)],
    'TNFSF12': [('TNFRSF12A',)],
}
RECEPTOR_SOURCES = {
    'scope': 'Manually curated conservative signalling routes for the original 14; other factors remain unreviewed, not receptor-negative. Gene subunits only; stoichiometry not inferred.',
    'BMP': 'https://reactome.org/content/detail/R-HSA-201451',
    'Activin': 'https://reactome.org/content/detail/R-HSA-1181153',
    'IFNG': 'https://reactome.org/content/detail/R-HSA-629623',
    'IL4_IL13': 'https://reactome.org/content/detail/R-HSA-449670',
    'TGFB': 'https://reactome.org/content/detail/R-HSA-170834',
    'IL22': 'https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId=1903',
    'TRAIL': 'https://reactome.org/content/detail/R-HSA-75238',
    'TWEAK': 'https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=5067',
    'INHBA_note': 'INHBA is the activin A homodimer subunit; expression does not establish processing or active ligand.',
    'TGFB_note': 'Latent TGF beta activation is not established by TGFB transcript detection.',
    'IL4_note': 'Only the epithelial type II IL4R/IL13RA1 route is used here; type I IL4R/IL2RG is not an epithelial substitute.'
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def significant(padj, fc):
    return math.isfinite(padj) and 0 <= padj < .05 and math.isfinite(fc) and abs(fc) > 1


def parse_workbook(source, factors, panel):
    """Parse once; retain signed DE records, diagnose every skipped row and duplicate sheet."""
    audits, records = [], []
    with ZipFile(source) as z:
        strings = string_table(z)
        for name, xmlpath in sheets(z):
            header = None
            counts = Counter()
            perts = Counter()
            digest = hashlib.sha256()
            columns = None
            for index, row in enumerate(cells_by_row(z, xmlpath, strings)):
                # Hash values, not styles, XML formatting, worksheet IDs or cell references.
                digest.update(json.dumps(sorted(row.items()), separators=(',', ':')).encode())
                digest.update(b'\n')
                if header is None:
                    h = find_header(row)
                    if all(h[k] is not None for k in ('gene','padj','logfc','pert')):
                        header, columns = h, [row[k] for k in sorted(row)]
                    elif index >= 25:
                        raise ValueError(f'{name}: required gene/padj/logFC/pert header not found')
                    continue
                counts['data_rows'] += 1
                factor = factor_name(row.get(header['pert'], ''))
                perts[factor] += 1
                if factor == 'NOCYTOMIX':
                    counts['control_rows'] += 1
                    continue
                if factor not in factors:
                    counts['unknown_factor_rows'] += 1
                    continue
                gene = str(row.get(header['gene'], '')).strip().upper()
                if not gene or gene == 'NA':
                    counts['missing_gene_rows'] += 1
                    continue
                counts['accepted_factor_gene_rows'] += 1
                padj = fnumber(row.get(header['padj'], ''))
                fc = fnumber(row.get(header['logfc'], ''))
                if not math.isfinite(padj) or not 0 <= padj <= 1:
                    counts['invalid_or_missing_padj_rows'] += 1
                    continue
                if padj < .05:
                    counts['fdr_significant_rows'] += 1
                    records.append({'factor':factor, 'sheet':name, 'gene':gene,
                                    'padj':padj, 'log2fc':fc, 'in_panel':int(gene in panel),
                                    'strong':int(significant(padj, fc))})
                    if significant(padj, fc):
                        counts['strong_rows'] += 1
                if not math.isfinite(fc):
                    counts['missing_log2fc_rows'] += 1
            if header is None:
                raise ValueError(f'{name}: no header')
            a = {'sheet':name, 'xml_path':xmlpath, 'columns':';'.join(columns),
                 'value_sha256':digest.hexdigest(), 'n_distinct_treatments':len(perts)}
            for key in ('data_rows','control_rows','unknown_factor_rows','missing_gene_rows',
                        'accepted_factor_gene_rows','invalid_or_missing_padj_rows',
                        'fdr_significant_rows','strong_rows','missing_log2fc_rows'):
                a[key] = counts[key]
            a['treatments'] = ';'.join(sorted(perts))
            audits.append(a)
            print(f'{name}: {counts["data_rows"]:,} rows, {counts["strong_rows"]:,} strong', flush=True)
    hashcounts = Counter(a['value_sha256'] for a in audits)
    for a in audits:
        a['duplicate_group'] = ';'.join(b['sheet'] for b in audits
            if b['value_sha256'] == a['value_sha256']) if hashcounts[a['value_sha256']] > 1 else ''
        a['use_for_conservative_candidates'] = int(not a['duplicate_group'])
    allowed = {a['sheet'] for a in audits if a['use_for_conservative_candidates']}
    if not records or not allowed:
        raise ValueError('No usable experimental response evidence')
    return audits, records, allowed


def factor_tables(factors, records, panel, allowed):
    rows, bycell = [], []
    for factor, info in sorted(factors.items()):
        allrec = [r for r in records if r['factor'] == factor]
        rec = [r for r in allrec if r['sheet'] in allowed]
        fdr = {r['gene'] for r in rec}
        strong = {r['gene'] for r in rec if r['strong']}
        covered = strong & panel
        ncell = 0
        for sheet in sorted({r['sheet'] for r in allrec}):
            sr = [r for r in allrec if r['sheet'] == sheet and r['strong']]
            genes = {r['gene'] for r in sr}
            overlap = genes & panel
            ncell += int(sheet in allowed and len(overlap) >= 5)
            bycell.append({'factor':factor, 'receiver_sheet':sheet,
                'source_usable':int(sheet in allowed), 'n_fdr_genes':len({r['gene'] for r in allrec if r['sheet']==sheet}),
                'n_strong_genes':len(genes), 'n_strong_genes_in_panel':len(overlap),
                'n_up_in_panel':len({r['gene'] for r in sr if r['log2fc']>1 and r['in_panel']}),
                'n_down_in_panel':len({r['gene'] for r in sr if r['log2fc']< -1 and r['in_panel']}),
                'coverage':len(overlap)/len(genes) if genes else 0})
        complex_ = CANONICAL_RECEPTORS.get(factor, [])
        complete = [r for r in complex_ if set(r) <= panel]
        row = {'factor':factor, 'ligand_gene':info['gene_symbol'],
               'ligand_in_panel':int(info['gene_symbol'] in panel),
               'catalogue_joint_lr_eligible':int(info['full_panel_lr_eligible']),
               'canonical_receptor_reviewed':int(factor in CANONICAL_RECEPTORS),
               'canonical_complexes_in_panel':';'.join('_'.join(r) for r in complete),
               'n_fdr_genes':len(fdr), 'n_strong_genes':len(strong),
               'n_strong_genes_in_panel':len(covered), 'response_coverage':len(covered)/len(strong) if strong else 0,
               'n_usable_receiver_types_with_5_targets':ncell,
               'all_sheets_n_strong_in_panel':len({r['gene'] for r in allrec if r['strong'] and r['in_panel']}),
               'conservative_panel_candidate':int(bool(info['gene_symbol'] in panel and complete and ncell))}
        rows.append(row)
    return rows, bycell


def complete_positive(positive, complexes):
    import numpy as np
    result = np.zeros(len(next(iter(positive.values()))), dtype=bool)
    for subunits in complexes:
        if all(g in positive for g in subunits):
            result |= np.logical_and.reduce([positive[g] for g in subunits])
    return result


def scan_expression(root, rows, records, allowed, min_total_umis, chunksize):
    """Aggregate only: never export cells. Targets are detection evidence, not activation scores."""
    import numpy as np
    import pandas as pd
    from specimens import SPECIMENS
    selected = [r for r in rows if r['conservative_panel_candidate']]
    programs = {(r['factor'], r['sheet']): sorted({s['gene'] for s in records
        if s['factor']==r['factor'] and s['sheet']==r['sheet'] and s['strong'] and s['in_panel']})
        for r in records if r['sheet'] in allowed and r['strong'] and r['in_panel']
        and any(x['factor']==r['factor'] for x in selected)}
    program_sets = {k:set(v) for k,v in programs.items()}
    gene_stats, specimen_stats, complex_stats, program_stats = [], [], [], []
    qc = Counter()
    for prefix, specs in SPECIMENS.items():
        path = root/'counts'/(prefix+'_exprMat_file.csv.gz')
        fov_to_geo = {f:geo for geo,cond,lo,hi in specs for f in range(lo,hi+1)}
        gene_counts, ligand_counts, rec_counts, comp_counts = {}, Counter(), Counter(), Counter()
        prog_counts, joint_counts = Counter(), Counter()
        print(f'Scanning {path.name}', flush=True)
        for chunk in pd.read_csv(path, chunksize=chunksize):
            if list(chunk.columns[:3]) != ['fov','cell_ID','cell']:
                raise ValueError('Unexpected expression header')
            mask = chunk['fov'].isin(fov_to_geo)
            chunk = chunk.loc[mask]
            if chunk.empty:
                continue
            numeric = chunk.iloc[:,3:]
            # Exclude technical controls from QC if present.
            numeric = numeric.loc[:,[c for c in numeric if not technical_target(c)]]
            passed = numeric.sum(axis=1) >= min_total_umis
            chunk = chunk.loc[passed]
            if chunk.empty:
                continue
            data = chunk.iloc[:,3:].to_numpy() > 0
            ix = {g.upper():i for i,g in enumerate(chunk.columns[3:])}
            positive = {g:data[:,i] for g,i in ix.items()}
            geos = chunk['fov'].map(fov_to_geo).to_numpy()
            for geo,_,_,_ in specs:
                selector = geos == geo
                qc[geo] += int(selector.sum())
                counts = data[selector].sum(axis=0).astype(int)
                gene_counts[geo] = gene_counts.get(geo, np.zeros(data.shape[1], dtype=np.int64)) + counts
            for row in selected:
                factor = row['factor']
                complexes = [tuple(x.split('_')) for x in row['canonical_complexes_in_panel'].split(';')]
                receiver = complete_positive(positive, complexes)
                ligand = positive[row['ligand_gene']]
                for geo,_,_,_ in specs:
                    selector = geos==geo
                    ligand_counts[geo,factor] += int(ligand[selector].sum())
                    rec_counts[geo,factor] += int(receiver[selector].sum())
                for complex_ in complexes:
                    cp = complete_positive(positive, [complex_])
                    for geo,_,_,_ in specs:
                        comp_counts[geo,factor,'_'.join(complex_)] += int(cp[geos==geo].sum())
                for key, genes in programs.items():
                    if key[0] != factor:
                        continue
                    detection = np.column_stack([positive[g] for g in genes]).sum(axis=1) >= 5
                    for geo,_,_,_ in specs:
                        selector = geos==geo
                        prog_counts[geo,key] += int(detection[selector].sum())
                        joint_counts[geo,key] += int((detection & receiver)[selector].sum())
        for geo,condition,_,_ in specs:
            if not qc[geo]:
                raise ValueError(f'No QC cells in {geo}')
            for gene,i in sorted(ix.items()):
                if any(gene in genes for genes in program_sets.values()) or any(gene==r['ligand_gene'] for r in selected):
                    gene_stats.append({'geo':geo,'gene':gene,'n_qc_cells':qc[geo],
                        'positive_cells':int(gene_counts[geo][i]),'detection_fraction':float(gene_counts[geo][i]/qc[geo])})
            for row in selected:
                f = row['factor']
                specimen_stats.append({'geo':geo,'condition':condition,'factor':f,'n_qc_cells':qc[geo],
                    'ligand_positive_cells':ligand_counts[geo,f], 'canonical_receptor_positive_cells':rec_counts[geo,f]})
                for c in row['canonical_complexes_in_panel'].split(';'):
                    complex_stats.append({'geo':geo,'factor':f,'complex':c,'n_qc_cells':qc[geo],
                                          'all_subunits_positive_cells':comp_counts[geo,f,c]})
                for key,genes in programs.items():
                    if key[0]!=f:continue
                    detected = sum(gene_counts[geo][ix[g]]>0 for g in genes)
                    program_stats.append({'geo':geo,'factor':f,'source_receiver_sheet':key[1],
                        'n_panel_targets':len(genes), 'n_targets_detected_in_specimen':detected,
                        'cells_with_5_targets_detected':prog_counts[geo,key],
                        'canonical_receptor_and_5_targets_same_cell':joint_counts[geo,key],
                        'receiver_cell_identity':'unannotated_tissue_cells'})
        print(f'QC: {dict(qc)}',flush=True)
    return gene_stats, specimen_stats, complex_stats, program_stats, dict(qc)


def save_table(out, name, rows, headers=None):
    if headers is None and not rows:
        raise ValueError(f'{name}: unexpected empty result; schema must be explicit')
    fields = headers or list(rows[0])
    writecsv(out/name, rows, fields)
    with (out/name).open(newline='') as f:
        rd = csv.DictReader(f)
        reread = list(rd)
        if rd.fieldnames != fields or len(reread) != len(rows) or any(None in r for r in reread):
            raise ValueError(f'{name}: CSV round-trip failed')
    return {'file':name, 'rows':len(rows), 'columns':len(fields),
            'bytes':(out/name).stat().st_size, 'sha256':sha256(out/name)}


def audit(args):
    out = args.out.expanduser()
    out.mkdir(parents=True, exist_ok=True)
    if (out/'response_overlap_summary.json').exists():
        raise FileExistsError('Use a fresh output directory to preserve previous runs')
    panel, sourcefiles = read_panel(args.root)
    prior = read_factor_metadata(args.previous)
    factors = {factor_name(r['factor_reagent_name']):r for r in prior}
    if len(factors) != len(prior) or len(factors) != 79:
        raise ValueError('Expected 79 unique canonical factors')
    audits, records, allowed = parse_workbook(args.supplement, factors, panel)
    rows, bycell = factor_tables(factors, records, panel, allowed)
    evidence = [dict(r, direction='up' if r['log2fc']>0 else 'down',
                     source_usable=int(r['sheet'] in allowed)) for r in records if r['strong']]
    tables = []
    for name, content in [('supplement_sheet_audit.csv', audits),
                           ('response_overlap_by_factor.csv',rows),
                           ('response_overlap_by_celltype.csv',bycell),
                           ('signed_response_genes.csv',evidence)]:
        tables.append(save_table(out,name,content))
    candidates = [dict(r) for r in rows if r['conservative_panel_candidate']]
    qc = {}
    if not args.panel_only:
        gs, ss, cs, ps, qc = scan_expression(args.root,rows,records,allowed,args.min_total_umis,args.chunksize)
        for name,content in [('target_gene_detection_by_specimen.csv',gs),
                              ('candidate_expression_by_specimen.csv',ss),
                              ('canonical_receptor_detection_by_specimen.csv',cs),
                              ('response_detection_by_specimen.csv',ps)]:
            tables.append(save_table(out,name,content))
        for r in candidates:
            f = r['factor']
            r['ligand_positive_cells'] = sum(x['ligand_positive_cells'] for x in ss if x['factor']==f)
            r['canonical_receptor_positive_cells'] = sum(x['canonical_receptor_positive_cells'] for x in ss if x['factor']==f)
            r['specimens_with_ligand_receptor_and_5_targets'] = sum(
                any(x['geo']==geo and x['factor']==f and x['ligand_positive_cells']>0
                    and x['canonical_receptor_positive_cells']>0 for x in ss)
                and any(x['geo']==geo and x['factor']==f and x['canonical_receptor_and_5_targets_same_cell']>0 for x in ps)
                for geo in qc)
    tables.append(save_table(out,'provisional_compatibility_candidates.csv', candidates,
                            list(candidates[0]) if candidates else list(rows[0])))
    # Independent recomputation of counts from exported tables is in validate_audit.py.
    summary = {'panel_gene_count_excluding_controls':len(panel), 'factor_treatments_audited':len(rows),
        'factors_with_any_fdr_gene':sum(r['n_fdr_genes']>0 for r in rows),
        'factors_with_at_least_2_strong_genes':sum(r['n_strong_genes']>=2 for r in rows),
        'factors_with_5_panel_strong_genes':sum(r['n_strong_genes_in_panel']>=5 for r in rows),
        'conservative_candidates':len(candidates), 'sheets_inspected':len(audits),
        'sheets_used_conservatively':sorted(allowed),
        'duplicate_sheets_excluded':[a['sheet'] for a in audits if not a['use_for_conservative_candidates']],
        'qc_cells_by_specimen':qc, 'n_qc_cells':sum(qc.values()),
        'thresholds':{'padj_lt':.05,'abs_log2fc_gt':1,'minimum_targets':5,'minimum_total_counts':args.min_total_umis},
        'source_files':[{'path':str(p),'sha256':sha256(p)} for p in [args.supplement,args.previous,args.catalog,*map(Path,sourcefiles)]],
        'code_sha256':sha256(Path(__file__)),
        'specimen_mapping_sha256':sha256(Path(__file__).with_name('specimens.py')),
        'limitations':[
            'Both duplicated sheets are excluded; stem/goblet identity unresolved. Other sheets retain published labels.',
            'Strong responses use both signs; paper FDR-only responder counts and up-only DEG counts are different definitions.',
            'No transcript detected is not proof of biological absence; >0 counts are permissive detection.',
            'Target detection is not signed programme activation, dose calibration or validated signalling.',
            'Tissue receivers lack cell-type annotation; co-detection is not proof of an epithelial response.',
            'Canonical receptors manually reviewed only for the original 14; other factors remain unreviewed.',
            'Spatial adjacency, ligand processing, receptor protein, specificity and held-out validation remain outstanding.',
            'Organoid effects are contextual additions versus cytomix controls; not normal tissue baseline effects.'
        ], 'tables':tables}
    (out/'receptor_sources.json').write_text(json.dumps(RECEPTOR_SOURCES,indent=2)+'\n')
    (out/'response_overlap_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('tables','source_files','limitations')},indent=2),flush=True)


def main():
    p = argparse.ArgumentParser(description='Conservative, signed organoid response audit; existing inputs only')
    base = Path.home()/'Real_Data/spatial_ibd_sheaf'
    p.add_argument('--root',type=Path,default=base/'pilot_6k')
    p.add_argument('--previous',type=Path,default=base/'outputs/organoid_79_audit'/FACTOR_FILE)
    p.add_argument('--supplement',type=Path,default=base/'metadata/Capeling_Supplementary_Data_7.xlsx')
    p.add_argument('--catalog',type=Path,default=base/'metadata/liana_omni_resource.csv')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--panel-only',action='store_true')
    p.add_argument('--no-download',action='store_true',help='Accepted for old CLI compatibility; this audit never downloads')
    p.add_argument('--min-total-umis',type=int,default=50)
    p.add_argument('--chunksize',type=int,default=512)
    audit(p.parse_args())

if __name__ == '__main__':
    main()
