"""Private receipt/hash/array QC after transfer. No fitting or observed-mask claims."""
from pathlib import Path
import argparse,json,hashlib,datetime,sys,tomllib
import numpy as np,pandas as pd,h5py
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE));from reconcile_native_identities import primary_identities

def main():
 ap=argparse.ArgumentParser();ap.add_argument('receipt',type=Path);a=ap.parse_args()
 base=ROOT/'data/worm_behavior_sheaf_readiness_v1';dest=base/'recordings'
 receipt=json.loads(a.receipt.read_text(encoding='utf-8-sig'))
 if receipt.get('algorithm','').upper().replace('-','')!='SHA256':raise ValueError('SHA256 receipt required')
 rows=receipt['files'];names=[r['filename'] for r in rows]
 allow=set(pd.read_csv(HERE/'REMAINING_FILES.csv')['file'])
 if len(names)!=len(set(names)) or not set(names)<=allow:raise ValueError('Unexpected or duplicate transfer filenames')
 for r in rows:
  p=dest/r['filename']
  if p.is_symlink() or not p.is_file():raise ValueError('Missing or symlink recording')
  if p.stat().st_size!=r['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256'].lower():raise ValueError('Transfer integrity failure')
 native=base/'sample_ingestion_v2/native_identity';nr=tomllib.loads((native/'native_receipt_PRIVATE.toml').read_text());labels=native/'labels_PRIVATE.tsv'
 if hashlib.sha256(labels.read_bytes()).hexdigest()!=nr['labels_sha256']:raise ValueError('Native TSV hash mismatch')
 if hashlib.sha256((dest/'NeuroPAL_labels.jld2').read_bytes()).hexdigest()!=nr['source_sha256']:raise ValueError('Native source hash mismatch')
 identities=primary_identities(pd.read_csv(labels,sep='\t',keep_default_na=False));results=[];coverage={}
 cfg=json.loads((HERE.parent/'READINESS_PROTOCOL.json').read_text())
 sys.path.insert(0,str(HERE.parent));from qc_primitives import valid_windows
 for p in sorted(dest.glob('*-cudata.h5')):
  animal=p.name.removesuffix('-cudata.h5');ids=identities[identities.animal==animal]
  with h5py.File(p,'r') as f:
   t=f['timing/timestamp_confocal'][:];x=f['gcamp/trace_array_original'][:]
   if x.ndim!=2 or x.shape[0]!=len(t) or not np.isfinite(t).all() or not (np.diff(t)>0).all():raise ValueError('Invalid temporal axes')
   neurons=[]
   for _,identity in ids.iterrows():
    i=int(identity.trace_row_1based)-1
    if not 0<=i<x.shape[1]:raise ValueError('Native row exceeds actual trace bounds')
    signal=x[:,i];fraction=float(np.isfinite(signal).mean());windows=int(valid_windows(signal,t,cfg['history_samples_primary']).sum())
    eligible=bool(fraction>=cfg['minimum_finite_trace_fraction'] and np.nanstd(signal)>0 and windows>=cfg['minimum_valid_short_history_windows_per_eligible_neuron'])
    neurons.append({'label':identity.label,'class':identity.neuron_class,'finite_fraction':fraction,'numerical_windows':windows,'eligible_processed_state':eligible})
   behaviour={}
   for key in ['velocity','reversal_vec','head_angle','pumping','angular_velocity','dist_from_cuCenter']:
    b=f['behavior/'+key][:]
    if len(b)!=len(t):raise ValueError('Behaviour length mismatch')
    behaviour[key]={'finite_fraction':float(np.isfinite(b).mean()),'std':float(np.nanstd(b))}
   coverage[animal]=sorted({r['class'] for r in neurons if r['eligible_processed_state']})
   results.append({'animal':animal,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'ROIs':x.shape[1],'frames':len(t),'neural_gap_count':int((np.diff(t)>cfg['timestamp_gap_multiplier']*np.median(np.diff(t))).sum()),'identities':neurons,'behaviour':behaviour,'genuine_observation_mask_available':False})
 stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out=base/'export_recovery_v3'/('cohort_audit_'+stamp+'_PRIVATE.json')
 out.write_text(json.dumps({'receipt_files_verified':len(rows),'animals':results,'class_coverage':coverage,'interpretation':'Numerical processed-state QC only; no original observation masks'},indent=2)+'\n')
 print(json.dumps({'receipt_files_verified':len(rows),'actual_recordings_inspected':len(results),'genuine_observation_masks_available':False,'private_results_saved':True,'full_readiness_not_automatically_declared':True}))
if __name__=='__main__':main()
