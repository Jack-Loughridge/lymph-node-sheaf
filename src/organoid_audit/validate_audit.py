#!/usr/bin/env python3
"""Independently check all CSV manifests, summary counts and gene unions.
Use --workbook for a second XLSX implementation (openpyxl read-only).
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import math
import re
import audit_organoid_response_genes as a


def read(path):
    with path.open(newline='') as f:
        rd=csv.DictReader(f)
        rows=list(rd)
        assert rd.fieldnames and all(None not in r for r in rows), path
        return rows


def validate(out, workbook=False):
    summary=json.loads((out/'response_overlap_summary.json').read_text())
    for m in summary['tables']:
        path=out/m['file']
        assert a.sha256(path)==m['sha256'], path
        rows=read(path)
        assert len(rows)==m['rows'], path
        assert len(rows[0])==m['columns'] if rows else True
        assert path.stat().st_size==m['bytes'], path
    rows=read(out/'response_overlap_by_factor.csv')
    byfactor={r['factor']:r for r in rows}
    assert len(byfactor)==len(rows)==summary['factor_treatments_audited']==79
    assert summary['factors_with_any_fdr_gene']==sum(int(r['n_fdr_genes'])>0 for r in rows)
    assert summary['factors_with_at_least_2_strong_genes']==sum(int(r['n_strong_genes'])>=2 for r in rows)
    assert summary['factors_with_5_panel_strong_genes']==sum(int(r['n_strong_genes_in_panel'])>=5 for r in rows)
    candidates=read(out/'provisional_compatibility_candidates.csv')
    assert len(candidates)==summary['conservative_candidates']
    assert {r['factor'] for r in candidates}=={r['factor'] for r in rows if int(r['conservative_panel_candidate'])}
    for r in candidates:
        assert int(r['ligand_in_panel']) and r['canonical_complexes_in_panel']
        assert int(r['n_usable_receiver_types_with_5_targets'])>0
    evidence=read(out/'signed_response_genes.csv')
    unique=defaultdict(set)
    for r in evidence:
        assert float(r['padj'])<.05 and abs(float(r['log2fc']))>1
        assert r['direction']==('up' if float(r['log2fc'])>0 else 'down')
        if int(r['source_usable']) and int(r['in_panel']):unique[r['factor']].add(r['gene'])
    for f,r in byfactor.items():assert len(unique[f])==int(r['n_strong_genes_in_panel']), f
    bycell=read(out/'response_overlap_by_celltype.csv')
    for f,r in byfactor.items():
        assert int(r['n_usable_receiver_types_with_5_targets'])==sum(
            x['factor']==f and int(x['source_usable']) and int(x['n_strong_genes_in_panel'])>=5 for x in bycell)
    ss=read(out/'candidate_expression_by_specimen.csv')
    cs=read(out/'canonical_receptor_detection_by_specimen.csv')
    ps=read(out/'response_detection_by_specimen.csv')
    gs=read(out/'target_gene_detection_by_specimen.csv')
    assert len(ss)==4*len(candidates)
    assert sum(summary['qc_cells_by_specimen'].values())==summary['n_qc_cells']
    for r in ss:
        n=summary['qc_cells_by_specimen'][r['geo']]
        assert int(r['n_qc_cells'])==n
        for key in ('ligand_positive_cells','canonical_receptor_positive_cells'):assert 0<=int(r[key])<=n
        comps=[x for x in cs if (x['geo'],x['factor'])==(r['geo'],r['factor'])]
        assert max(int(x['all_subunits_positive_cells']) for x in comps)<=int(r['canonical_receptor_positive_cells'])
        assert int(r['canonical_receptor_positive_cells'])<=sum(int(x['all_subunits_positive_cells']) for x in comps)
    for r in ps:
        ssr=next(x for x in ss if (x['geo'],x['factor'])==(r['geo'],r['factor']))
        assert 0<=int(r['canonical_receptor_and_5_targets_same_cell'])<=min(int(ssr['canonical_receptor_positive_cells']),int(r['cells_with_5_targets_detected']))
        assert 0<=int(r['n_targets_detected_in_specimen'])<=int(r['n_panel_targets'])
    detected={(x['geo'],x['gene']):int(x['positive_cells']) for x in gs}
    for x in ps:
        targets={r['gene'] for r in evidence if r['factor']==x['factor']
                 and r['sheet']==x['source_receiver_sheet'] and int(r['in_panel']) and int(r['source_usable'])}
        assert len(targets)==int(x['n_panel_targets'])
        assert sum(detected[x['geo'],g]>0 for g in targets)==int(x['n_targets_detected_in_specimen'])
    assert len({(x['geo'],x['gene']) for x in gs})==len(gs)
    for r in gs:
        assert 0<=int(r['positive_cells'])<=int(r['n_qc_cells'])
        assert abs(float(r['detection_fraction'])-int(r['positive_cells'])/int(r['n_qc_cells']))<1e-12
    for r in candidates:
        f=r['factor']
        for field in ('ligand_positive_cells','canonical_receptor_positive_cells'):
            assert int(r[field])==sum(int(x[field]) for x in ss if x['factor']==f)
        assert int(r['specimens_with_ligand_receptor_and_5_targets'])==sum(
            any(x['geo']==geo and x['factor']==f and int(x['ligand_positive_cells'])>0 and int(x['canonical_receptor_positive_cells'])>0 for x in ss)
            and any(x['geo']==geo and x['factor']==f and int(x['canonical_receptor_and_5_targets_same_cell'])>0 for x in ps)
            for geo in summary['qc_cells_by_specimen'])
    report={'all_table_checks':'passed','table_count':len(summary['tables']),
            'summary_count_checks':'passed','candidates':len(candidates)}
    if workbook:
        import openpyxl
        source=Path(summary['source_files'][0]['path'])
        wb=openpyxl.load_workbook(source, read_only=True, data_only=True)
        expected=read(out/'supplement_sheet_audit.csv')
        old_path=Path(summary['source_files'][1]['path']).parent.parent/'organoid_response_audit/response_overlap_by_factor.csv'
        old={a.factor_name(r['factor']):r for r in read(old_path)}
        panel,_=a.read_panel(Path(summary['source_files'][3]['path']).parent.parent)
        allsig, allstrong, consigs, constr = (defaultdict(set) for _ in range(4))
        sheets=[]
        for ws in wb:
            counter=Counter()
            iterator=ws.iter_rows(values_only=True)
            cols=next(iterator)
            assert cols==('gene','baseMean','log2FoldChange','lfcSE','pvalue','padj','pert'), cols
            for gene,base,fc,se,p,padj,pert in iterator:
                counter['data_rows']+=1
                factor=a.factor_name(pert)
                if factor not in byfactor:continue
                # Published workbook stores many scientific-notation values as TEXT.
                try:padj=float(padj)
                except (ValueError,TypeError):padj=float('nan')
                try:fc=float(fc)
                except (ValueError,TypeError):fc=float('nan')
                if math.isfinite(padj) and 0<=padj<.05:
                    allsig[factor].add(str(gene).strip().upper())
                    counter['fdr_significant_rows']+=1
                    if ws.title in summary['sheets_used_conservatively']:consigs[factor].add(str(gene).strip().upper())
                    if math.isfinite(fc) and abs(fc)>1:
                        counter['strong_rows']+=1
                        allstrong[factor].add(str(gene).strip().upper())
                        if ws.title in summary['sheets_used_conservatively']:constr[factor].add(str(gene).strip().upper())
            original=next(x for x in expected if x['sheet']==ws.title)
            for key in ('data_rows','fdr_significant_rows','strong_rows'):assert counter[key]==int(original[key]), (ws.title,key,counter[key],original[key])
            sheets.append({'sheet':ws.title,**counter})
            print('Independent openpyxl check passed:',ws.title,flush=True)
        wb.close()
        for f,r in byfactor.items():
            assert len(consigs[f])==int(r['n_fdr_genes']), f
            assert len(constr[f])==int(r['n_strong_genes']), f
            assert len(allstrong[f]&panel)==int(old[f]['response_genes_in_panel']), f
            legacy_sig={g for g in allsig[f] if re.fullmatch(r'[A-Z0-9.-]{2,35}',g)}
            assert len(legacy_sig)==int(old[f]['significant_fdr_genes']), (f,len(legacy_sig),old[f]['significant_fdr_genes'])
        report['independent_workbook_parser']='passed'
        report['original_factor_counts_reproduced']='79/79 using the original gene-name filter'
        report['genes_rejected_by_legacy_filter']={f:sorted(g for g in gs if not re.fullmatch(r'[A-Z0-9.-]{2,35}',g)) for f,gs in allsig.items() if any(not re.fullmatch(r'[A-Z0-9.-]{2,35}',g) for g in gs)}
        report['independent_sheet_counts']=sheets
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('out',type=Path)
    p.add_argument('--workbook',action='store_true')
    args=p.parse_args()
    validate(args.out,args.workbook)
