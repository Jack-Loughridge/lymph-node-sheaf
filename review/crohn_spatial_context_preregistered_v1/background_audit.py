"""Audit controls without changing any prior criterion, result or exposure call."""
from pathlib import Path
import json,csv,hashlib
import numpy as np,h5py
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
SOURCE=ROOT/'data/crohn_xenium/scp3259/xenium5k_scVI_CellCharter_batch2_reannot_SCP.h5ad'
OLD=ROOT/'review/crohn_xenium_biological_validation_v1'
def per_control_background(total_probe_counts,n_cells,n_controls):
 if n_cells<=0 or n_controls<=0:raise ValueError('Verified positive denominators required')
 return total_probe_counts/(n_cells*n_controls)
def main():
 before={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in OLD.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
 before.update({name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['review/signalling_sheaf/BIOLOGICAL_ACCEPTANCE_CRITERIA.md','src/signalling_sheaf/MODEL_SPECIFICATION.txt']})
 snapshot=HERE/'FROZEN_ORIGINAL_SHA256.json'
 if not snapshot.exists(): snapshot.write_text(json.dumps(before,indent=2)+'\n')
 frozen=json.loads(snapshot.read_text())
 for p,d in frozen.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==d,p
 z=np.load(ROOT/'data/crohn_xenium/biological_validation_v1/cell_metadata_arrays_PRIVATE.npz');qc=z['qc']
 with h5py.File(SOURCE,'r') as f:
  sc=f['obs/biosample_id/codes'][:];probes=f['obs/control_probe_counts'][:].astype(np.float64);codewords=f['obs/control_codeword_counts'][:].astype(np.float64)
  genes=f['var/_index'].asstr()[:];controls=[g for g in genes if 'negcontrol' in g.lower() or 'negativecontrol' in g.lower()]
  rows=[]
  for k in range(16):
   for cohort,mask in [('all_exported_cells',sc==k),('primary_QC_cells',(sc==k)&qc)]:
    n=int(mask.sum());probe_sum=float(probes[mask].sum());code_sum=float(codewords[mask].sum());combined=probes[mask]+codewords[mask]
    # Conditional denominator 20 from general 10x guidance, NOT verified for this run.
    mean=per_control_background(probe_sum,n,20)
    rows.append(dict(sample_alias=f'S{k+1:02d}',cohort=cohort,cells=n,negative_probe_sum=int(probe_sum),negative_codeword_sum=int(code_sum),mean_probe_total_per_cell=probe_sum/n,mean_codeword_total_per_cell=code_sum/n,assumed_control_count=20,control_count_verified=False,conditional_mean_per_control_per_cell=mean,conditional_false_positive_per_exported_5001_genes=mean*len(genes),old_pooled_99_9_percentile=float(np.quantile(combined,.999)),old_minimum_count_floor=float(max(2,1+np.quantile(combined,.999))),interpretation='Conditional aggregate QC estimate; not gene/cell positivity threshold'))
 with (HERE/'background_metrics.csv').open('w') as out:
  w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 summary=dict(official_equivalent=False,negative_control_features_retained_in_gene_matrix=len(controls),panel_control_denominator_verified=False,probe_and_codeword_counts_available_separately=True,per_probe_profiles_available=False,original_files_verified_unchanged=len(frozen),run_specific_panel_json_and_analysis_summary_missing=True)
 (HERE/'BACKGROUND_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
 for p,d in frozen.items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==d,p
 print(json.dumps(summary))
if __name__=='__main__':main()
