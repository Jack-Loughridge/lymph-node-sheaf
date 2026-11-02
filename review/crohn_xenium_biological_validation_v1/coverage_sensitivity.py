"""Prespecified raw-floor sensitivity plus diagnostic author-label populations.
No response genes, patient grouping, predictive fitting or receptor imputation.
The original structural marker screen is unavailable on this panel; do not treat
its empty receiver set as evidence that fibroblasts are absent.
"""
import csv,json
from pathlib import Path
import numpy as np,h5py
from scipy.spatial import cKDTree
import spatial_audit as a


def main():
 meta=np.load(a.LOCAL/'cell_metadata_arrays_PRIVATE.npz');mc=json.loads((a.LOCAL/'selected_gene_columns_PRIVATE.json').read_text());v=np.memmap(a.LOCAL/'selected_raw_counts.int32.memmap',dtype='int32',mode='r',shape=tuple(mc['shape']));gi={g:i for i,g in enumerate(mc['genes'])}
 with h5py.File(a.SOURCE) as f:
  types=f['obs/annot_new/categories'].asstr()[:];broad=f['obs/megagut_lvl1_majority_voting/categories'].asstr()[:]
  aliases=[f'S{i+1:02d}' for i in range(16)];cross=json.loads((a.LOCAL/'sample_patient_crosswalk_PRIVATE.json').read_text())
  lc=meta['annotation_codes'];bc=meta['broad_codes'];sc=meta['sample_codes'];xy=meta['coords'];area=meta['cell_area'];conf=meta['confidence'];neg=meta['negative_control_counts'];tot=meta['raw_totals'];ng=meta['raw_genes'];qc=meta['qc']
  mescode=int(np.where(broad=='Mesenchymal')[0][0]);fib=a.label_mask(np.array(['fibroblast' in x.lower() for x in types]),lc)&(bc==mescode)
  immune=np.isin(bc,[i for i,b in enumerate(broad) if b in ['B and B plasma','Myeloid','T and NK cells']])
  inf=a.label_mask(np.array([x=='M COL7A1hi inflammatory fibroblast' for x in types]),lc)
  rows=[];marker=[];gene_stats=[];qc_counts=[]
  for k in range(16):
   sample=sc==k;qmask=sample&qc;D=float(np.median(2*np.sqrt(area[qmask]/np.pi)));floor=max(2,int(np.quantile(neg[qmask],.999,method='higher'))+1)
   populations={'author_fibroblast_labels':qmask&fib,'author_fibroblast_labels_conf_ge_08':qmask&fib&(conf>=.8),'author_inflammatory_fibroblast_label':qmask&fib&inf,'all_author_mesenchymal_labels':qmask&(bc==mescode)}
   for mincounts in [20,50,100]:
    qm=sample&(tot>=mincounts)&(ng>=20)&np.isfinite(xy).all(1)&(area>0)
    qc_counts.append(dict(sample_alias=aliases[k],minimum_counts=mincounts,QC_cells=int(qm.sum()),fibroblast_labelled=int((qm&fib).sum())))
   for pop,mask in populations.items():
    mm=v[mask];coords=xy[mask];ctrl=neg[mask];N=len(mm)
    for gene in ['FAP','PDPN','PDGFRA','COL4A1','COL7A1','DPT','PI16']:
     if gene not in gi:continue
     marker.append(dict(sample_alias=aliases[k],population=pop,gene=gene,cells=N,detected_ge1=int((mm[:,gi[gene]]>=1).sum())))
    for programme,route in a.ROUTES.items():
     cols=[gi[g] for g in route['receptor']];lig=gi[route['ligand']]
     raw_comp={thr:np.all(mm[:,cols]>=thr,axis=1) for thr in [1,2,3]}
     comp=np.all(mm[:,cols]>=floor,axis=1)&np.all(mm[:,cols]>ctrl[:,None],axis=1)
     send=qmask&immune&(v[:,lig]>=floor)&(v[:,lig]>neg)
     distances=cKDTree(xy[send]).query(coords,k=1)[0] if send.any() and N else np.full(N,np.inf)
     for radius_mult in [1.5,2.5,4.]:
      exposed=comp&(distances<=radius_mult*D)
      rows.append(dict(sample_alias=aliases[k],condition_from_name=cross[k]['status_from_explicit_sample_name'],population=pop,programme=programme,radius_multiplier=radius_mult,eligible_cells=N,raw_receptor_ge1=int(raw_comp[1].sum()),raw_receptor_ge1_fraction=float(raw_comp[1].mean()) if N else None,raw_receptor_ge2=int(raw_comp[2].sum()),raw_receptor_ge3=int(raw_comp[3].sum()),screened_receptor_competent=int(comp.sum()),immune_senders=int(send.sum()),exposed=int(exposed.sum()),exposed_fraction=float(exposed.mean()) if N else None,control_floor=floor))
    # Additional empirical detection/quantiles, not external response-map evidence.
    for gene,col in gi.items():
     vals=mm[:,col];nz=vals[vals>0]
     gene_stats.append(dict(sample_alias=aliases[k],population=pop,gene=gene,cells=N,detected_ge1=int((vals>0).sum()),mean_counts=float(vals.mean()) if N else None,p99_counts=float(np.quantile(vals,.99)) if N else None,max_counts=int(vals.max()) if N else None))
  for name,data in [('diagnostic_receiver_coverage.csv',rows),('panel_marker_concordance.csv',marker),('QC_sensitivity.csv',qc_counts),('selected_gene_distributions.csv',gene_stats)]:
   with (a.HERE/name).open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
  print('Diagnostic coverage saved; original marker-screen gate remains NOT EVALUABLE. No patient analyses.')
if __name__=='__main__':main()
