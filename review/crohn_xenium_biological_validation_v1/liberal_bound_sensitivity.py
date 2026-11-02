"""No-background >=1-count sensitivity is a liberal plausibility bound, not exposure call.
Compare prespecified QC floors and radii; no patient aggregation or outcome fit.
"""
import csv,json
import numpy as np
from scipy.spatial import cKDTree
import h5py
import spatial_audit as a

def main():
 d=np.load(a.LOCAL/'cell_metadata_arrays_PRIVATE.npz');g=json.loads((a.LOCAL/'selected_gene_columns_PRIVATE.json').read_text());x=np.memmap(a.LOCAL/'selected_raw_counts.int32.memmap',dtype='int32',mode='r',shape=tuple(g['shape']));gi={k:i for i,k in enumerate(g['genes'])}
 with h5py.File(a.SOURCE) as f:
  labels=f['obs/annot_new/categories'].asstr()[:];broad=f['obs/megagut_lvl1_majority_voting/categories'].asstr()[:]
 lc=d['annotation_codes'];bc=d['broad_codes'];sc=d['sample_codes'];coords=d['coords'];area=d['cell_area'];conf=d['confidence'];tot=d['raw_totals'];ng=d['raw_genes'];mes=int(np.where(broad=='Mesenchymal')[0][0])
 fib_refined=a.label_mask(np.array(['fibroblast' in s.lower() for s in labels]),lc);fib=fib_refined&(bc==mes);imm=np.isin(bc,[i for i,b in enumerate(broad) if b in ['B and B plasma','Myeloid','T and NK cells']]);out=[]
 for k in range(16):
  for mincounts in [20,50,100,200]:
   qc=(sc==k)&(tot>=mincounts)&(ng>=20)&np.isfinite(coords).all(1)&(area>0);diam=float(np.median(2*np.sqrt(area[qc]/np.pi)))
   for pop,rec in [('author_fibroblast_labels',qc&fib),('author_fibroblast_labels_conf_ge_08',qc&fib&(conf>=.8)),('all_refined_fibroblast_labels_without_broad_consensus',qc&fib_refined),('all_refined_fibroblast_labels_conf_ge08_without_broad_consensus',qc&fib_refined&(conf>=.8))]:
    rx=x[rec];rxy=coords[rec]
    for programme,route in a.ROUTES.items():
     rc=[gi[g] for g in route['receptor']];li=gi[route['ligand']]
     for threshold in [1,2,3]:
      comp=np.all(rx[:,rc]>=threshold,axis=1);send=qc&imm&(x[:,li]>=threshold)
      dist=cKDTree(coords[send]).query(rxy,k=1)[0] if send.any() and len(rxy) else np.full(len(rxy),np.inf)
      for multiple in [1.5,2.5,4.]:
       exposed=comp&(dist<=diam*multiple)
       out.append(dict(sample_alias=f'S{k+1:02d}',minimum_panel_counts=mincounts,population=pop,programme=programme,detection_count=threshold,radius_multiplier=multiple,receivers=len(rx),competent=int(comp.sum()),competent_fraction=float(comp.mean()) if len(comp) else None,immune_ligand_senders=int(send.sum()),spatial_candidate_receivers=int(exposed.sum()),candidate_fraction=float(exposed.mean()) if len(exposed) else None,interpretation='No control correction; liberal sensitivity only'))
 with (a.HERE/'liberal_detection_radius_QC_sensitivity.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
 print('Liberal detection/QC/radius diagnostics:',len(out),'aggregate rows')
if __name__=='__main__':main()
