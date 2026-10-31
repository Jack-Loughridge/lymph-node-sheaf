"""Read-only source inspection, strict composite-key joins and tiny local cache.

Cache contains selected cell observations and belongs under ignored data/, not Git.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd
from allocation import INPUT_GENES

MARKERS={
 'epithelial':['EPCAM','KRT8','KRT18','KRT19','CLDN4'],
 'immune':['PTPRC','CD3D','CD3E','LST1','TYROBP'],
 'stromal':['COL1A1','COL1A2','DCN','PDGFRA'],
 'cycling':['MKI67','PCNA','TOP2A','CCNA2','MCM5'],
 'absorptive':['CA2','SLC26A2','GUCA2B','KRT20','FABP1'],
 'goblet':['MUC2','TFF3'],
}
ANNOTATION_GENES=sorted({g for v in MARKERS.values() for g in v})
SELECTION=[
 ('Cos1400R355RRNA6K4u1716',37,'GSM9345416'),
 ('Cos1400R355RRNA6K4u1716',25,'GSM9345417'),
 ('Cos1380R387RRNA6K4u2716',21,'GSM9345418'),
 ('Cos1380R387RRNA6K4u2716',6,'GSM9345419'),
]
SHEET_TYPES={'A_TA_DEGs':'TA_like','D_Colonocyte_DEGs':'colonocyte_like'}
FACTOR_PROGRAMME={'IFNG':'IFNG','TNF':'TNF','BMP2':'BMP','BMP6':'BMP'}


def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()


def strict_join(expr,meta):
 """Slide+FOV+cell, never numeric cell_ID alone; no multiplicative/partial join."""
 keys=['prefix','fov','cell']
 if expr[keys].isna().any().any() or meta[keys].isna().any().any():raise ValueError('Missing join keys')
 if expr.duplicated(keys).any() or meta.duplicated(keys).any():raise ValueError('Duplicate cell keys')
 joined=expr.merge(meta,on=keys,how='left',validate='one_to_one',indicator=True,suffixes=('','_meta'))
 if not (joined['_merge']=='both').all():raise ValueError('Unmatched expression-to-metadata cell keys')
 if 'cell_ID_meta' in joined and not (joined['cell_ID'].to_numpy()==joined['cell_ID_meta'].to_numpy()).all():
  raise ValueError('Numeric cell IDs disagree across join')
 return joined.drop(columns='_merge')


def labels(counts,genes):
 ix={g:i for i,g in enumerate(genes)}
 missing=set(ANNOTATION_GENES)-set(ix)
 if missing:raise ValueError('Missing annotation markers: '+str(sorted(missing)))
 detected={k:(counts[:,[ix[g] for g in group]]>0).sum(axis=1) for k,group in MARKERS.items()}
 epi=(detected['epithelial']>=2)&(detected['immune']<2)&(detected['stromal']<2)
 result=np.full(len(counts),'other_or_mixed',dtype='U30')
 result[epi]='epithelial_unresolved'
 ta=epi&(detected['cycling']>=2)&(detected['absorptive']<2)&(detected['goblet']<2)
 colon=epi&(detected['absorptive']>=2)&(detected['cycling']<2)&(detected['goblet']<2)
 result[ta]='TA_like';result[colon]='colonocyte_like'
 return result,detected


def signatures(path):
 rows=[]; per= {}
 with Path(path).open(newline='') as f:
  for r in csv.DictReader(f):
   if r['factor'] not in FACTOR_PROGRAMME:continue
   r=dict(r)
   reason=[]
   if r['source_usable']!='1':reason.append('ambiguous_duplicate_sheet')
   if r['in_panel']!='1':reason.append('not_on_panel')
   if r['gene'] in INPUT_GENES:reason.append('ligand_or_receptor_input')
   if r['gene'] in ANNOTATION_GENES:reason.append('annotation_marker')
   r['exclusion_reason']=';'.join(reason)
   r['receiver_proxy']=SHEET_TYPES.get(r['sheet'],'unmatched_source')
   rows.append(r)
   if not reason:
    key=(r['receiver_proxy'],r['factor'])
    group=per.setdefault(key,{})
    if r['gene'] in group:raise ValueError('Repeated source gene/factor/type')
    group[r['gene']]=float(np.sign(float(r['log2fc'])))
 weights={}
 for c in SHEET_TYPES.values():
  for p in ('IFNG','TNF','BMP'):
   factors=['BMP2','BMP6'] if p=='BMP' else [p]
   w={}
   for factor in factors:
    g=per.get((c,factor),{})
    if not g:raise ValueError('Missing signed programme '+str((c,factor)))
    for gene,sign in g.items():w[gene]=w.get(gene,0)+sign/(len(g)*len(factors))
   weights[c,p]={g:v for g,v in w.items() if abs(v)>1e-14}
 return rows,weights


def prepare(root,signature,cache):
 root=Path(root);cache=Path(cache); signature=Path(signature)
 rows,weights=signatures(signature)
 wanted=sorted(set(INPUT_GENES)|set(ANNOTATION_GENES)|{g for w in weights.values() for g in w})
 sources=[signature]
 for prefix in dict.fromkeys(x[0] for x in SELECTION):
  sources += [root/'counts'/(prefix+'_exprMat_file.csv.gz'),root/'metadata'/(prefix+'_metadata_file.csv')]
 fingerprint=[{'path':str(p.resolve()),'sha256':digest(p)} for p in sources]
 info_path=cache.with_suffix('.json')
 if cache.exists():
  if not info_path.exists():raise ValueError('Cache missing provenance')
  info=json.loads(info_path.read_text())
  if info['sources']!=fingerprint or info['genes']!=wanted:raise ValueError('Stale cache; use a new cache path')
  print('Verified existing selected-FOV cache',flush=True)
  return np.load(cache,allow_pickle=False),info
 arrays=[];inspect=[];joins=[]
 for prefix in dict.fromkeys(x[0] for x in SELECTION):
  matrix=root/'counts'/(prefix+'_exprMat_file.csv.gz')
  meta_path=root/'metadata'/(prefix+'_metadata_file.csv')
  meta=pd.read_csv(meta_path)
  inspect.append({'file':str(meta_path),'rows':len(meta),'columns':list(meta.columns),
   'published_annotation_columns':[c for c in meta if any(k in c.lower() for k in ('celltype','cell_type','annotation','cluster','predicted_label'))],
   'duplicate_composite_keys':int(meta.duplicated(['fov','cell']).sum()),
   'interpretation':'No validated epithelial subtype column; fluorescence channels are not subtype labels'})
  chosen={f:geo for pre,f,geo in SELECTION if pre==prefix}
  meta=meta.loc[meta.fov.isin(chosen)].copy();meta['prefix']=prefix
  import gzip
  with gzip.open(matrix,'rt') as stream:panel=next(csv.reader(stream))[3:]
  missing=set(wanted)-set(panel)
  if missing:raise ValueError('Required features absent from panel: '+str(sorted(missing)))
  inspect[-1]['marker_panel_audit']={
    'initial_absorptive_candidates_missing':sorted(set(['CA1','CA4','SLC26A3','GUCA2A','AQP8','MS4A12','ALPI','FABP2','CEACAM7'])-set(panel)),
    'initial_goblet_candidates_missing':sorted(set(['SPINK4','CLCA1','ZG16'])-set(panel)),
    'used_marker_sets':MARKERS,
    'warning':'Available absorptive markers are not specific enough for validated colonocyte labels; use proxies only'}
  expr=[];selected_before=0;selected_after=0
  print('Reading selected FOVs from '+matrix.name,flush=True)
  for chunk in pd.read_csv(matrix,chunksize=512):
   if list(chunk.columns[:3])!=['fov','cell_ID','cell']:raise ValueError('Unexpected expression header')
   chunk=chunk.loc[chunk.fov.isin(chosen)].copy()
   if chunk.empty:continue
   biological=[g for g in chunk.columns[3:] if not re.fullmatch(r'(Negative|SystemControl)[0-9]+',g)]
   if len(biological)!=6175:raise ValueError('Unexpected biological panel')
   total=chunk[biological].sum(axis=1)
   selected_before+=len(chunk)
   chunk=chunk.loc[total>=50].copy();chunk['total']=total.loc[chunk.index]
   selected_after+=len(chunk)
   chunk['prefix']=prefix
   expr.append(chunk[['prefix','fov','cell','cell_ID','total',*wanted]])
  if not expr:raise ValueError('Empty selected FOV inputs')
  expr=pd.concat(expr,ignore_index=True)
  joined=strict_join(expr,meta)
  if joined[['CenterX_global_px','CenterY_global_px','Area']].isna().any().any():raise ValueError('Missing spatial values')
  for f,geo in chosen.items():
   frame=joined.loc[joined.fov==f].sort_values('cell').copy()
   if not len(frame):raise ValueError('No QC cells in FOV')
   arrays.append(dict(counts=frame[wanted].to_numpy(float), genes=np.array(wanted),
    total=frame['total'].to_numpy(float),fov=frame.fov.to_numpy(int),
    geo=np.full(len(frame),geo),prefix=np.full(len(frame),prefix),cell=frame.cell.to_numpy(str),
    coords=frame[['CenterX_global_px','CenterY_global_px']].to_numpy(float),
    areas=frame['Area'].to_numpy(float),
    panck=frame['Mean.PanCK'].to_numpy(float),cd45=frame['Mean.CD45'].to_numpy(float),cd68=frame['Mean.CD68'].to_numpy(float)))
   joins.append({'geo':geo,'prefix':prefix,'fov':f,'metadata_cells':int((meta.fov==f).sum()),
                 'qc_expression_cells':len(frame),'joined_cells':len(frame),'unmatched_expression_cells':0})
  print(f'Joined {selected_after}/{selected_before} selected cells after biological QC',flush=True)
 result={k:(arrays[0][k] if k=='genes' else np.concatenate([x[k] for x in arrays],axis=0)) for k in arrays[0]}
 labs,_=labels(result['counts'],wanted);result['labels']=labs
 cache.parent.mkdir(parents=True,exist_ok=True)
 np.savez_compressed(cache,**result)
 info={'sources':fingerprint,'selection':SELECTION,'selection_reason':'Previously chosen median metadata-cell-count FOV per specimen; fixed before response analysis',
       'genes':wanted,'inspections':inspect,'joins':joins,'cache_contains_cell_observations':'Local ignored data/ only; do not commit or include in review archive'}
 info_path.write_text(json.dumps(info,indent=2)+'\n')
 return np.load(cache,allow_pickle=False),info

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--root',type=Path,default=Path.home()/'Real_Data/spatial_ibd_sheaf/pilot_6k')
 p.add_argument('--signatures',type=Path,default=Path('outputs/organoid_response_audit_corrected/signed_response_genes.csv'))
 p.add_argument('--cache',type=Path,default=Path('data/signalling_sheaf_inputs.npz'))
 args=p.parse_args();prepare(args.root,args.signatures,args.cache)
