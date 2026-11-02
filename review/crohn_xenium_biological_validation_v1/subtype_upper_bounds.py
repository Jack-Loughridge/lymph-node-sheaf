"""Exhaustive author-refined fibroblast-subtype receptor bounds, no subtype selection.
These ignore ligand/background/distance and are necessary, not sufficient gates.
"""
import json,csv
import numpy as np,h5py
import spatial_audit as a

def main():
 d=np.load(a.LOCAL/'cell_metadata_arrays_PRIVATE.npz');g=json.loads((a.LOCAL/'selected_gene_columns_PRIVATE.json').read_text());x=np.memmap(a.LOCAL/'selected_raw_counts.int32.memmap',mode='r',dtype='int32',shape=tuple(g['shape']));gi={k:i for i,k in enumerate(g['genes'])}
 with h5py.File(a.SOURCE) as f:labels=f['obs/annot_new/categories'].asstr()[:]
 out=[]
 for label_code,label in enumerate(labels):
  if 'fibroblast' not in label.lower():continue
  for minimum in [20,100,200]:
   eligible=(d['annotation_codes']==label_code)&(d['raw_totals']>=minimum)&(d['raw_genes']>=20)&(d['confidence']>=.8)&np.isfinite(d['coords']).all(1)&(d['cell_area']>0)
   for programme,route in a.ROUTES.items():
    nr=[];co=[]
    for k in range(16):
     r=x[eligible&(d['sample_codes']==k)];nr.append(len(r));co.append(int(np.all(r[:,[gi[g] for g in route['receptor']]]>=1,axis=1).sum()))
    fractions=[c/n if n else 0 for c,n in zip(co,nr)]
    out.append(dict(author_subtype=label,minimum_panel_counts=minimum,programme=programme,total_eligible=sum(nr),total_complete_receptor_ge1=sum(co),specimens_ge20percent=sum(v>=.2 for v in fractions),max_fraction=max(fractions),necessary_upper_bound_gate=bool(sum(co)>=100 and sum(v>=.2 for v in fractions)>=3)))
 with (a.HERE/'all_fibroblast_subtype_receptor_bounds.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
 print('Necessary receptor-only gates by programme:',{p:sum(r['necessary_upper_bound_gate'] for r in out if r['programme']==p) for p in a.ROUTES})
if __name__=='__main__':main()
