"""Additional PANEL-COMPATIBLE annotation diagnostics after marker-availability check.
Not a replacement eligibility rule. Aggregate all-gene fibroblast detection in a
second sparse streaming pass; no clinical mapping or response fitting.
"""
import json,csv
from collections import Counter
import numpy as np,h5py
import spatial_audit as a


def expected_lineage(label):
 s=label.lower()
 if any(x in s for x in ['macrophage','monocyte','mast']) or 'dc' in s:return 'Myeloid'
 if 'fibroblast' in s or 'muscle' in s or 'pericyte' in s:return 'Mesenchymal'
 if 'endotheli' in s or s.startswith('ec '):return 'Endothelial'
 if 'enterocyte' in s or 'goblet' in s or 'epitheli' in s:return 'Epithelial'
 if 'neuron' in s or 'glial' in s:return 'Neural'
 if s.startswith('b_'):return 'B and B plasma'
 if s.startswith('t_') or s.startswith('th_'):return 'T and NK cells'
 return None


def main():
 z=np.load(a.LOCAL/'cell_metadata_arrays_PRIVATE.npz');sc=z['sample_codes'];lc=z['annotation_codes'];bc=z['broad_codes'];qc=z['qc'];conf=z['confidence']
 with h5py.File(a.SOURCE) as f:
  labels=f['obs/annot_new/categories'].asstr()[:];broad=f['obs/megagut_lvl1_majority_voting/categories'].asstr()[:];genes=f['var/_index'].asstr()[:].tolist();ng=len(genes);n=len(sc)
  predicted=np.full(len(lc),None,dtype=object);predicted[lc>=0]=np.asarray([expected_lineage(s) for s in labels],object)[lc[lc>=0]];observed=broad[bc];matches=predicted==observed
  fib=a.label_mask(np.array(['fibroblast' in s.lower() for s in labels]),lc)&(observed=='Mesenchymal');inflammatory=a.label_mask(np.array([s=='M COL7A1hi inflammatory fibroblast' for s in labels]),lc)
  pops=[qc&fib,qc&fib&(conf>=.8),qc&fib&inflammatory];pnames=['author_fibroblast_labels','author_fibroblast_labels_conf_ge_08','author_inflammatory_fibroblast_label']
  extra_genes=[g for g in ['FAP','PDPN'] if g in genes];extra=np.memmap(a.LOCAL/'extra_marker_counts.int32.memmap',mode='w+',dtype='int32',shape=(n,len(extra_genes)))
  det=np.zeros((3,16,ng),np.int64);det2=np.zeros_like(det);sums=np.zeros((3,16,ng),np.float64)
  for start in range(0,n,25000):
   stop=min(n,start+25000);m=a.read_csr_rows(f['layers/counts'],start,stop,ng);m.sum_duplicates();m.eliminate_zeros();m.sort_indices()
   extra[start:stop]=m[:,[genes.index(g) for g in extra_genes]].toarray().astype(np.int32)
   for k in np.unique(sc[start:stop]):
    for j,pop in enumerate(pops):
     mm=m[(sc[start:stop]==k)&pop[start:stop]]
     det[j,k]+=np.asarray((mm>0).sum(0)).ravel();det2[j,k]+=np.asarray((mm>=2).sum(0)).ravel();sums[j,k]+=np.asarray(mm.sum(0)).ravel()
  extra.flush();np.savez_compressed(a.LOCAL/'all_gene_fibroblast_detection_PRIVATE.npz',genes=np.asarray(genes),population=np.asarray(pnames),detected_ge1=det,detected_ge2=det2,counts_sum=sums)
  rows=[];norm=[]
  for k in range(16):
   qq=qc&(sc==k);ff=qq&fib
   row=dict(sample_alias=f'S{k+1:02d}',QC_cells=int(qq.sum()),refined_broad_concordant=int(matches[qq].sum()),refined_annotation_coverage=float((lc[qq]>=0).mean()),refined_annotation_cells=int((qq&(lc>=0)).sum()),refined_broad_concordance=float(matches[qq&(lc>=0)].mean()),fibroblast_cells=int(ff.sum()))
   for j,g in enumerate(extra_genes):row[g+'_positive']=int((extra[ff,j]>0).sum());row[g+'_positive_fraction']=float((extra[ff,j]>0).mean())
   row['FAP_or_PDPN_positive']=int(np.any(extra[ff]>0,axis=1).sum());row['FAP_and_PDPN_positive']=int(np.all(extra[ff]>0,axis=1).sum())
   rows.append(row)
   start=int(np.flatnonzero(sc==k)[0]);stop=min(start+1000,n);raw=a.read_csr_rows(f['layers/counts'],start,stop,ng);analysis=a.read_csr_rows(f['X'],start,stop,ng)
   raw.sort_indices();analysis.sort_indices();same=np.array_equal(raw.indptr,analysis.indptr) and np.array_equal(raw.indices,analysis.indices)
   denom=f['obs/total_counts'][start:stop];expect=np.log1p(raw.data*10000/np.repeat(denom,np.diff(raw.indptr)))
   norm.append(dict(sample_alias=f'S{k+1:02d}',rows_checked=stop-start,aligned_sparse_structure=same,CP10k_log1p_allclose=bool(same and np.allclose(analysis.data,expect,rtol=1e-5,atol=1e-5))))
  with (a.HERE/'annotation_concordance.csv').open('w') as out:
   w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  with (a.HERE/'normalisation_checks.csv').open('w') as out:
   w=csv.DictWriter(out,fieldnames=list(norm[0]));w.writeheader();w.writerows(norm)
  summary=dict(initial_structural_markers=a.FIBRO_MARKERS,initial_structural_markers_measured=[g for g in a.FIBRO_MARKERS if g in genes],initial_ge2_marker_eligibility_status='NOT EVALUABLE: only one prescribed marker is measured; empty mask is not absence of fibroblasts',diagnostic_panel_markers=extra_genes,independent_orthogonal_validation=False,broad_refined_field_concordance_among_annotated=float(matches[qc&(lc>=0)].mean()),refined_annotation_coverage_in_QC=float((lc[qc]>=0).mean()),normalised_X_checked_rows=sum(x['rows_checked'] for x in norm),normalisation_checks_all_pass=all(x['CP10k_log1p_allclose'] for x in norm))
  a.save(a.HERE/'ANNOTATION_NORMALISATION_SUMMARY.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
