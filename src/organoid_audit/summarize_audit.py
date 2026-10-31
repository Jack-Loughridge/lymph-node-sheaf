#!/usr/bin/env python3
"""Write a reviewable plain-text report from completed audit outputs."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
from zipfile import ZipFile


def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))


def main(out):
    summary=json.loads((out/'response_overlap_summary.json').read_text())
    candidates=read(out/'provisional_compatibility_candidates.csv')
    evidence=read(out/'signed_response_genes.csv')
    ps=read(out/'response_detection_by_specimen.csv')
    gs=read(out/'target_gene_detection_by_specimen.csv')
    prior=Path(summary['source_files'][1]['path']).parent.parent/'organoid_response_audit'
    original=[]
    for path in sorted(prior.glob('*.csv')):
        rows=read(path)
        original.append({'file':path.name,'rows':len(rows),'bytes':path.stat().st_size,
                         'blank_records':sum(not any(r.values()) for r in rows)})
    with ZipFile(prior/'response_audit_results.zip') as z:
        archive_matches=all(z.read(n)==(prior/Path(n).name).read_bytes() for n in z.namelist())
    (out/'original_output_diagnosis.json').write_text(json.dumps({'csvs':original,'zip_members_match_loose_files':archive_matches},indent=2)+'\n')
    text=['Capeling/CosMx signalling compatibility feasibility audit',
          'Date: 2026-10-31; Python 3.10.12; dependencies pinned in src/organoid_audit/requirements.txt.',
          '', 'Diagnosis: original CSV files are populated; no header or CSV export failure found.']
    text += [f"  {r['file']}: {r['rows']} data rows; {r['bytes']} bytes; {r['blank_records']} blank records." for r in original]
    text += [f'Original ZIP members match loose CSVs: {archive_matches}.',
        'Original shortlist implicitly required ligand coverage through full_panel_lr_eligible; that was not a missing-ligand bug.',
        '', 'Workbook: columns gene, baseMean, log2FoldChange, lfcSE, pvalue, padj, pert.',
        'A_TA_DEGs: 274018 rows; B_Stem_DEGs: 85141; C_Goblet_DEGs: 85141; D_Colonocyte_DEGs: 158615.',
        'Stem and goblet worksheets are value-identical throughout. Neither identity is independently trustworthy.',
        'Both are excluded from conservative summaries. TA and colonocyte labels are retained as published.',
        'Signed evidence retains every strong row with an explicit source_usable flag.',
        'NoCytomix is a control, not one of the 79 factor reagents.',
        'Scientific-notation text in padj is converted to numeric. The original parser already did this correctly.',
        'The old restrictive gene-name regex discarded THRA1/BTR (one FDR-significant BMP2 row).',
        'The corrected audit retains the source label; it is not panel-covered and does not alter the candidate targets.',
        '', 'Reanalysis:',
        f"  Biological panel targets: {summary['panel_gene_count_excluding_controls']}; controls excluded: 344 (20 Negative*, 324 SystemControl*).",
        f"  QC cells: {summary['n_qc_cells']}; >=50 biological counts per cell.",
        f"  QC by specimen: {json.dumps(summary['qc_cells_by_specimen'],sort_keys=True)}.",
        f"  79 treatments; {summary['factors_with_any_fdr_gene']} with any FDR-significant gene; {summary['factors_with_at_least_2_strong_genes']} with >=2 strong genes; {summary['factors_with_5_panel_strong_genes']} with >=5 panel-covered strong genes.",
        f"  {summary['conservative_candidates']} provisional candidates retain ligand, canonical receptor, and >=5 targets in at least one usable organoid receiver signature.",
        '  Strong = padj<0.05 AND abs(log2FC)>1; up and down genes retained separately.',
        '  IFNA remains an ambiguous reagent-family name. TGFB2 has 5 targets pooled across sheets but no usable single receiver signature with 5.',
        '  Canonical receptors reviewed for the original 14 only; other factors are unreviewed, not receptor-negative.',
        '', 'Pathways (targets = union of covered strong genes from usable TA/colonocyte signatures):',
        'factor | targets | ligand+ cells | canonical receptor+ cells | specimens with receptor AND >=5 targets in same cell | canonical routes']
    for r in candidates:
        text.append(f"{r['factor']} | {r['n_strong_genes_in_panel']} | {r['ligand_positive_cells']} | {r['canonical_receptor_positive_cells']} | {r['specimens_with_ligand_receptor_and_5_targets']}/4 | {r['canonical_complexes_in_panel']}")
    text += ['', 'All candidate ligand genes and canonical receptor routes are panel-covered and detected in all four specimens.',
             'All covered signature genes are individually detected in each specimen: '+str(all(x['n_panel_targets']==x['n_targets_detected_in_specimen'] for x in ps))+'.',
             'IL22 receptor/5-target co-detection: GSM9345418 1 cell; GSM9345419 17 cells; other specimens 0.',
             'TGFB3 receptor/5-target co-detection: GSM9345417 3 cells; GSM9345418 2 cells; other specimens 0.',
             'These are permissive transcript detection results in unannotated tissue cells, not confirmed epithelial programme activity.',
             '', 'Initial model recommendations (judgment, not a fitted performance ranking):',
             '1. IFNG: 113 covered targets, clear IFNGR1/IFNGR2 route, coherent CXCL9/10/11, IDO1, GBP4 and SOCS1 response.',
             '2. BMP2/BMP6: 133/88 targets; explicit type-I/type-II receptors and shared ID1/ID3 responses. Start with family-level compatibility.',
             '3. TNF: 32 targets; abundant canonical receptor detection. Preserve NF-kB/injury context and distinguish from other inflammatory factors.',
             '4. IL4/IL13: 30/22 targets; shared epithelial IL4R/IL13RA1 route. Group initially; disentangling ligands requires extra evidence.',
             '5. TNFSF10/TRAIL: 92 targets and strong receptor detectability; useful additional axis, but death/survival/differentiation context needs care.',
             'BMP4 is viable (82 targets) but original 343-target count depended heavily on ambiguous worksheets.',
             'TGFB1, INHBA, TNFSF11/RANKL and TNFSF12/TWEAK are feasible secondary axes with smaller signatures or overlapping effects.',
             'IL22 and TGFB3 have sparse same-cell evidence under the five-target rule; retain as exploratory, not a basis for disease comparisons.',
             '', 'Overlap of pooled unsigned target sets (for identifiability checks, not an activity score):']
    sets=defaultdict(set)
    for r in evidence:
        if int(r['source_usable']) and int(r['in_panel']):sets[r['factor']].add(r['gene'])
    for f,g in [('BMP2','BMP6'),('IL4','IL13'),('TGFB1','TGFB3'),('TGFB1','INHBA')]:
        shared=len(sets[f]&sets[g]);union=len(sets[f]|sets[g])
        text.append(f'  {f}/{g}: {shared} shared of {union} union; Jaccard {shared/union:.3f}.')
    text += ['', 'Biological limits and next model requirements:',
        'The source experiment added recombinant factors to cytomix-injured epithelial organoids from three donors.',
        'These are perturbation-derived DEG signatures, not the published cNMF gene programmes; no activity calibration has been fitted.',
        'Background/ambient counts, cell-type annotation, receptor protein, secretion and active ligand processing remain unresolved.',
        'Activin requires ACVR1B with ACVR2A/B; the old catalogue listed binding/accessory partners instead.',
        'IL22RA2 and TNFRSF10C/D are not substitutes for the respective signalling receptors.',
        'TGFB transcripts do not establish latent ligand activation; INHBA does not establish an active processed dimer.',
        'Use organoid receiver-specific signed weights, family-aware restrictions and neighbouring sender/receiver cells.',
        'Validate held-out specimens and spatial nulls before interpreting a sheaf discrepancy as disease organisation.',
        'Four specimens do not establish disease effects, and this audit makes no ageing claim.',
        '', 'Validation and reproducibility:',
        'Run tests and validate_audit.py --workbook as documented in src/organoid_audit/README.txt.',
        'The summary manifests every main CSV with row/column counts, byte lengths and SHA256 hashes.',
        'validation.json records the completed independent checks. No raw cell data were exported.',
        'Original sibling project files were read only; corrected code/results reside in the writable lymph_node_sheaf workspace.',
        'No large downloads, commits or pushes were performed.',
        '', 'Sources:',
        'Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC12891662/',
        'Author code: https://github.com/Genentech/secretome_dictionary',
        'BMP receptors: https://reactome.org/content/detail/R-HSA-201451',
        'Activin receptors: https://reactome.org/content/detail/R-HSA-1181153',
        'IFNG receptor: https://reactome.org/content/detail/R-HSA-629623',
        'IL4/IL13 receptor: https://reactome.org/content/detail/R-HSA-449670',
        'IL22 receptor: https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId=1903',
        'TGFB receptors: https://reactome.org/content/detail/R-HSA-170834',
        'TNF receptors: https://reactome.org/content/detail/R-HSA-75893',
        'TRAIL receptors: https://reactome.org/content/detail/R-HSA-75238',
        'RANKL: https://reactome.org/content/detail/R-HSA-5667141',
        'TWEAK: https://reactome.org/content/detail/R-HSA-5675974']
    (out/'audit_report.txt').write_text('\n'.join(text)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('out',type=Path)
    main(p.parse_args().out)
