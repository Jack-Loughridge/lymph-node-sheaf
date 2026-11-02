"""Chunked actual-data coverage audit; no full dense matrix or response fit.
All cell arrays and exact patient/sample crosswalks stay under ignored data/.
"""
import json,csv,re,hashlib,os,time
from pathlib import Path
from collections import Counter
import h5py,numpy as np,pandas as pd
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
SOURCE=ROOT/'data/crohn_xenium/scp3259/xenium5k_scVI_CellCharter_batch2_reannot_SCP.h5ad'
LOCAL=ROOT/'data/crohn_xenium/biological_validation_v1';LOCAL.mkdir(exist_ok=True)
CFG=json.loads((HERE/'AUDIT_PROTOCOL.json').read_text());ROUTES=CFG['primary_programmes']
MARKERS={
 'Mesenchymal':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','ACTA2','TAGLN','MYH11'],
 'Epithelial':['EPCAM','KRT8','KRT18','KRT19','CLDN4'],
 'Endothelial':['PECAM1','VWF','EMCN','KDR','RAMP2'],
 'B and B plasma':['PTPRC','CD79A','MS4A1','CD37','CD74','JCHAIN','MZB1'],
 'Myeloid':['PTPRC','LST1','TYROBP','CTSS','FCER1G','AIF1'],
 'T and NK cells':['PTPRC','CD3D','CD3E','TRAC','NKG7','GNLY','CD247'],
 'Neural':['SOX10','S100B','PLP1','TUBB3','SNAP25','NCAM1']}
FIBRO_MARKERS=['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA']


def save(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def decode(v):return v.decode() if isinstance(v,bytes) else str(v)
def category(f,name):
 g=f['obs/'+name];return g['categories'].asstr()[:],g['codes'][:]
def label_mask(category_flags,codes):
 codes=np.asarray(codes)
 if np.any(codes < -1) or np.any(codes >= len(category_flags)):raise ValueError('Invalid categorical code')
 return (codes>=0)&np.asarray(category_flags)[np.maximum(codes,0)]
def read_csr_rows(g,start,stop,ng):
 ptr=g['indptr'][start:stop+1].astype(np.int64);lo,hi=int(ptr[0]),int(ptr[-1]);ptr-=lo
 return csr_matrix((g['data'][lo:hi],g['indices'][lo:hi],ptr),shape=(stop-start,ng))
def marker_availability(measured,markers,minimum=2):
 available=sorted(set(measured)&set(markers))
 return dict(available=available,evaluable=len(available)>=minimum)
def co_detect(values,columns,threshold):
 return np.all(values[:,columns]>=threshold,axis=1)
def spatial_exposed(receiver_xy,sender_xy,radius):
 if not len(sender_xy):return np.zeros(len(receiver_xy),bool)
 d=cKDTree(sender_xy).query(receiver_xy,k=1)[0]
 return d<=radius


def map_samples(cats,workbook=None):
 if workbook is None:
  payload=json.loads((LOCAL/'workbook_linkage_PRIVATE.json').read_text())
  assert payload['source_sha256']==hashlib.sha256((ROOT/'data/crohn_xenium/261023_Xenium_metadata_linkage.xlsx').read_bytes()).hexdigest()
  wb=pd.DataFrame(payload['rows'])
 else:wb=workbook.copy()
 wb=wb[wb['Table S1 patient ID'].str.fullmatch(r'\d+(?:\.0)?')].copy()
 norm=lambda s:re.sub('[^a-z0-9]','',str(s).lower())
 rows=[]
 for i,s in enumerate(cats):
  slide=s.split('_')[0];cand=wb[wb['Xenium_slide'].str.zfill(7)==slide]
  if len(cand)>1:cand=cand[cand['Biobank ID'].map(lambda v:norm(v) in norm(s))]
  if len(cand)!=1:raise ValueError('Ambiguous specimen/workbook join; do not group patients')
  state='adjacent' if re.search(r'(adjacent|_adj$|_nor$|NOR$)',s) else 'affected' if re.search(r'(affected|_inf$|INF$)',s) else None
  if state is None:raise ValueError('Unresolved sample status')
  row=cand.iloc[0]
  rows.append(dict(biosample=s,biosample_code=i,alias=f'S{i+1:02d}',SCP_category_order=i,
    slide_id=slide,patient_id=row['Table S1 patient ID'],biobank_id=row['Biobank ID'],
    status_from_explicit_sample_name=state,biobank_token_used=len(wb[wb['Xenium_slide'].str.zfill(7)==slide])>1,
    S1_verified=False,mapping_status='PROVISIONAL_AUTHOR_WORKBOOK_AND_BIOSAMPLE_NAMES; S1 PENDING'))
 save(LOCAL/'sample_patient_crosswalk_PRIVATE.json',rows)
 return rows


def main():
 began=time.monotonic()
 with h5py.File(SOURCE,'r') as f:
  genes=f['var/_index'].asstr()[:].tolist();ng=len(genes);n=f['obs/x_centroid'].shape[0];assert len(set(genes))==ng
  samples,sc=category(f,'biosample_id');broad,bc=category(f,'megagut_lvl1_majority_voting');labels,lc=category(f,'annot_new');statuscats,statuscodes=category(f,'inflammation_status')
  cross=map_samples(samples);assert len(cross)==16
  coords=f['obsm/spatial'][:];centroids=np.c_[f['obs/x_centroid'][:],f['obs/y_centroid'][:]]
  assert coords.shape==(n,2) and np.array_equal(coords,centroids)
  area=f['obs/cell_area'][:];conf=f['obs/xenium4_conf_score'][:];control=f['obs/control_probe_counts'][:]+f['obs/control_codeword_counts'][:]
  cellcats,cc=category(f,'cell_id');composite=sc.astype(np.int64)*len(cellcats)+cc
  duplicate_composite=n-len(np.unique(composite));assert duplicate_composite==0
  idx=f['obs/_index'].asstr()[:];index_duplicates=n-len(set(idx));del idx,composite
  selected=sorted(set([g for v in MARKERS.values() for g in v]+FIBRO_MARKERS+[r['ligand'] for r in ROUTES.values()]+[g for r in ROUTES.values() for g in r['receptor']]+['OSM','OSMR','IL6ST','LIFR','TNFRSF1B','COL7A1','ADAMDEC1','COL4A1','DPT','SOCS3','CA12','CCL19','GREM1','MCAM','PI16','TNC','PDGFRB','MMP1','MMP3','IL11']))
  selected=[g for g in selected if g in genes];gj={g:i for i,g in enumerate(selected)};cols=[genes.index(g) for g in selected]
  # Dense ONLY for a small declared gene subset, not 5001 genes. Stored privately.
  path=LOCAL/'selected_raw_counts.int32.memmap';sel=np.memmap(path,dtype='int32',mode='w+',shape=(n,len(selected)))
  totals=np.zeros(n,np.float64);nexpr=np.zeros(n,np.int32);qc=np.zeros(n,bool)
  detection=np.zeros((16,ng),np.int64);gene_sums=np.zeros((16,ng),np.float64)
  negative=noninteger=nonfinite=0;rowsum_meta_diff=0;storedtot=f['obs/total_counts'][:]
  for start in range(0,n,25000):
   stop=min(n,start+25000);m=read_csr_rows(f['layers/counts'],start,stop,ng);m.sum_duplicates();m.eliminate_zeros();m.sort_indices()
   negative+=int((m.data<0).sum());nonfinite+=int((~np.isfinite(m.data)).sum());noninteger+=int((m.data!=np.round(m.data)).sum())
   if negative or nonfinite or noninteger:raise ValueError('Count layer contains unexpected noninteger/negative/nonfinite values')
   totals[start:stop]=np.asarray(m.sum(axis=1)).ravel();nexpr[start:stop]=m.getnnz(axis=1)
   q=(totals[start:stop]>=100)&(nexpr[start:stop]>=20)&np.isfinite(coords[start:stop]).all(axis=1)&np.isfinite(area[start:stop])&(area[start:stop]>0)
   qc[start:stop]=q;sel[start:stop]=m[:,cols].toarray().astype(np.int32)
   rowsum_meta_diff+=int((totals[start:stop]!=storedtot[start:stop]).sum())
   for k in np.unique(sc[start:stop]):
    mm=m[(sc[start:stop]==k)&q]
    detection[k]+=np.asarray((mm>0).sum(axis=0)).ravel();gene_sums[k]+=np.asarray(mm.sum(axis=0)).ravel()
   if start//25000%20==0:print('streamed_rows',stop,flush=True)
  sel.flush();np.savez_compressed(LOCAL/'gene_detection_by_specimen_PRIVATE.npz',detection=detection,sums=gene_sums,genes=np.asarray(genes),sample_aliases=np.asarray([r['alias'] for r in cross]))
  save(LOCAL/'selected_gene_columns_PRIVATE.json',dict(genes=selected,shape=[n,len(selected)],dtype='int32',source_layer='layers/counts'))
  np.savez_compressed(LOCAL/'cell_metadata_arrays_PRIVATE.npz',sample_codes=sc,broad_codes=bc,annotation_codes=lc,qc=qc,coords=coords,cell_area=area,confidence=conf,negative_control_counts=control,raw_totals=totals,raw_genes=nexpr)
  fibro_label=np.array(['fibroblast' in lab.lower() for lab in labels]);bmes=int(np.where(broad=='Mesenchymal')[0][0]);immune_classes={i for i,b in enumerate(broad) if b in ['B and B plasma','Myeloid','T and NK cells']}
  marker_fibro=(sel[:,[gj[g] for g in FIBRO_MARKERS if g in gj]]>0).sum(axis=1)
  is_fibro=label_mask(fibro_label,lc)&(bc==bmes);primary_receiver=qc&is_fibro&(conf>=.8)&(marker_fibro>=2)
  immune=qc&np.isin(bc,list(immune_classes));allresults=[];sample_summary=[];ann=[];gene_rows=[]
  for k,row in enumerate(cross):
   in_sample=sc==k;qmask=in_sample&qc;rec=in_sample&primary_receiver;allfib=qmask&is_fibro;mes=qmask&(bc==bmes)
   floor=max(2,int(np.quantile(control[qmask],.999,method='higher'))+1)
   diameter=float(np.median(2*np.sqrt(area[qmask]/np.pi)));rr=sel[rec];neg=control[rec]
   named=labels[lc[qmask&(lc>=0)]];sub_mask=(lc==int(np.where(labels=='M COL7A1hi inflammatory fibroblast')[0][0]))
   sample_summary.append(dict(sample_alias=row['alias'],condition_from_name=row['status_from_explicit_sample_name'],cells=int(in_sample.sum()),QC_cells=int(qmask.sum()),QC_fraction=float(qmask.sum()/in_sample.sum()),fibroblast_QC_cells=int(allfib.sum()),primary_receivers=int(rec.sum()),inflammatory_fibroblast_receivers=int((rec&sub_mask).sum()),mesenchymal_QC_cells=int(mes.sum()),radius_unit='micrometre_by_unscaled_Xenium_centroid_provenance',median_diameter_um=diameter,control_floor_counts=floor,control_positive_fraction=float((control[qmask]>0).mean()),centroid_x_range_um=[float(coords[qmask,0].min()),float(coords[qmask,0].max())],centroid_y_range_um=[float(coords[qmask,1].min()),float(coords[qmask,1].max())],source_inflammation_status_missing=int((qmask&(statuscodes==0)).sum())))
   for b,lin in enumerate(broad):
    mask=qmask&(bc==b);gg=[gj[g] for g in MARKERS[lin] if g in gj]
    corroborated=(sel[mask][:,gg]>0).sum(axis=1)>=2
    ann.append(dict(sample_alias=row['alias'],lineage=lin,cells=int(mask.sum()),marker_genes_measured=len(gg),cells_ge2_lineage_markers=int(corroborated.sum()),marker_corrob_fraction=float(corroborated.mean()) if len(corroborated) else None,marker_screen_status='NOT_EVALUABLE_PANEL_INSUFFICIENT' if len(gg)<2 else 'DIAGNOSTIC_TRANSCRIPT_CONCORDANCE_ONLY'))
   for prog,route in ROUTES.items():
    rc=[gj[g] for g in route['receptor']];lig=gj[route['ligand']]
    upper=co_detect(rr,rc,1);competent=co_detect(rr,rc,floor)&np.all(rr[:,rc]>neg[:,None],axis=1)
    send=qmask&immune&(sel[:,lig]>=floor)&(sel[:,lig]>control)
    sendxy=coords[send];recxy=coords[rec]
    for mult in [1.5,2.5,4.]:
     exposed=competent&spatial_exposed(recxy,sendxy,diameter*mult)
     # Bounding rectangle screen is descriptive; no exact tissue mask exists in this object.
     allxy=coords[qmask];margin=np.minimum.reduce([recxy[:,0]-allxy[:,0].min(),allxy[:,0].max()-recxy[:,0],recxy[:,1]-allxy[:,1].min(),allxy[:,1].max()-recxy[:,1]])
     interior=margin>=diameter*mult
     allresults.append(dict(sample_alias=row['alias'],condition_from_name=row['status_from_explicit_sample_name'],programme=prog,radius_multiplier=mult,radius_um=diameter*mult,eligible_receivers=int(rec.sum()),receptor_ge1_upper_bound=int(upper.sum()),upper_bound_fraction=float(upper.mean()) if len(upper) else None,complete_receptor_background_screen=int(competent.sum()),immune_ligand_senders=int(send.sum()),exposed_receivers=int(exposed.sum()),exposed_fraction=float(exposed.mean()) if len(exposed) else None,interior_eligible=int(interior.sum()),interior_exposed=int((exposed&interior).sum())))
    for gene in [route['ligand'],*route['receptor']]:
     v=sel[allfib,gj[gene]];iv=sel[qmask&immune,gj[gene]]
     gene_rows.append(dict(sample_alias=row['alias'],programme=prog,gene=gene,fibroblast_cells=len(v),fibroblast_detected_ge1=int((v>=1).sum()),fibroblast_detected_ge2=int((v>=2).sum()),immune_cells=len(iv),immune_detected_ge1=int((iv>=1).sum()),immune_detected_ge2=int((iv>=2).sum())))
  for name,rows in [('specimen_QC.csv',sample_summary),('lineage_marker_audit.csv',ann),('exposure_coverage.csv',allresults),('input_gene_detection.csv',gene_rows)]:
   with (HERE/name).open('w') as out:
    w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
  summary=dict(cells=n,genes=ng,nonzeros=int(f['layers/counts/data'].shape[0]),count_values_negative=negative,count_values_noninteger=noninteger,count_values_nonfinite=nonfinite,obs_total_differs_from_exported_panel_raw_sum=rowsum_meta_diff,QC_cells=int(qc.sum()),QC_fraction=float(qc.mean()),selected_dense_genes=len(selected),selected_cache_bytes=path.stat().st_size,composite_sample_cell_id_duplicates=duplicate_composite,obs_index_duplicates=index_duplicates,valid_coordinates=int(np.isfinite(coords).all(1).sum()),raw_cell_id_categories=len(cellcats),biosamples=len(samples),source_donor_id_categories=len(f['obs/donor_id/categories']),provisional_workbook_patients=len(set(r['patient_id'] for r in cross)),S1_mapping_verified=False,patient_grouped_analyses_performed=False,primary_receivers=int(primary_receiver.sum()),primary_marker_screen_status='NOT_EVALUABLE_PANEL_INSUFFICIENT',zero_primary_mask_is_not_absence_of_fibroblasts=True,refined_annotation_missing_cells=int((lc<0).sum()),refined_annotation_missing_QC_cells=int((qc&(lc<0)).sum()),refined_annotation_categories=len(labels),broad_annotation_categories=len(broad),counts_integer_verified=True,OSMR_present='OSMR' in genes,primary_genes_all_measured=all(g in genes for v in ROUTES.values() for g in [v['ligand'],*v['receptor']]),source_inflammation_status_missing=int((statuscodes==0).sum()),wall_seconds=time.monotonic()-began)
  save(HERE/'SPATIAL_AUDIT_SUMMARY.json',summary);print(json.dumps(summary,indent=2),flush=True)
  # Save compact masks for deterministic secondary detection/QC sensitivities, not public cell observations.
  np.savez_compressed(LOCAL/'receiver_masks_PRIVATE.npz',allfib=is_fibro,primary=primary_receiver,immune=immune,marker_fibro=marker_fibro)
if __name__=='__main__':main()
