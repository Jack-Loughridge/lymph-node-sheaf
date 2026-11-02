"""Read small transferred HDF5 recordings; never fit or synthesize observations.
Per-animal/neuron outputs stay ignored. Aggregate readiness is fail-closed.
Requires verified author-derived label TSV for identity coverage, not JLD2 guessing.
"""
from pathlib import Path
import argparse,csv,json,hashlib
from collections import defaultdict
import numpy as np,h5py
from qc_primitives import identity_join,valid_windows,common_core_exists
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
LOCAL=ROOT/'data/worm_behavior_sheaf_readiness_v1';RECORDINGS=LOCAL/'recordings'
def save(path,o):path.write_text(json.dumps(o,indent=2,allow_nan=False)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--labels-tsv',type=Path);ap.add_argument('--provenance-json',type=Path);args=ap.parse_args()
 paths=sorted(RECORDINGS.glob('*-cudata.h5'));labels=[]
 if args.labels_tsv:labels=list(csv.DictReader(args.labels_tsv.open(),delimiter='\t'))
 provenance=json.loads(args.provenance_json.read_text()) if args.provenance_json else {}
 rows=[];nr=[];files=[];schemas={};classes=defaultdict(set);individuals=defaultdict(set)
 for p in paths:
  animal=p.name.removesuffix('-cudata.h5')
  if p.stat().st_size>20_000_000:raise ValueError('Unexpected processed size; inspect first')
  files.append(dict(file=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
  with h5py.File(p,'r') as f:
   schema={}
   def visit(name,obj):
    if isinstance(obj,h5py.Dataset):schema[name]=dict(shape=list(obj.shape),dtype=str(obj.dtype))
   f.visititems(visit);schemas[animal]=schema
   g=f['gcamp/trace_array_original']
   if g.ndim!=2 or g.size>20_000_000:raise ValueError('Unexpected trace shape')
   x=g[:];t=f['timing/timestamp_confocal'][:].reshape(-1)
   if x.shape[1]!=len(t):raise ValueError('Trace/time axis mismatch')
   # No claims about units, sampling Hz or reconstruction from names alone.
   dt=np.diff(t);time_ok=bool(np.isfinite(t).all() and len(dt) and (dt>0).all())
   review=provenance.get(animal,{})
   segments=review.get('verified_segment_starts_zero_based',[])
   ids=identity_join([r for r in labels if r['animal']==animal],len(x))
   byrow={r['row']:r for r in ids}
   observed=np.isfinite(x)
   if review.get('observed_mask_dataset'):
    key=review['observed_mask_dataset'];m=f[key][:]
    if m.dtype!=bool or m.shape!=x.shape:raise ValueError('Provenance mask axes/type mismatch')
    observed &= m
   for i,a in enumerate(x):
    fraction=float(observed[i].mean());std=float(np.nanstd(a)) if np.isfinite(a).any() else None
    windows=int(valid_windows(a,t,10,observed[i],segments).sum()) if time_ok else 0
    eligible=bool(fraction>=.95 and windows>=200 and std is not None and std>0)
    ident=byrow.get(i,{})
    nr.append(dict(animal=animal,trace_row_1based=i+1,finite_fraction=float(np.isfinite(a).mean()),observed_candidate_fraction=fraction,valid_windows=windows,constant=(std==0),negative_values=int((a<0).sum()),label=ident.get('label',''),neuron_class=ident.get('neuron_class',''),confidence=ident.get('confidence',None),LR=ident.get('LR',''),numerically_eligible=eligible,genuine_observation_provenance_verified=bool(review.get('observed_signal_provenance_verified',False))))
    if eligible and ident.get('confidence',0)>=4:
     if ident.get('neuron_class'):classes[animal].add(ident['neuron_class'])
     if ident.get('label'):individuals[animal].add(ident['label'])
   behavior={}
   for key in ['velocity','reversal_vec','reversal_events','head_angle','pumping','dist_from_cuCenter','circle_radius']:
    if 'behavior/'+key in f:
     b=f['behavior/'+key][:] if f['behavior/'+key].shape else np.array(f['behavior/'+key][()])
     behavior[key]=dict(shape=list(b.shape),finite_fraction=float(np.isfinite(b).mean()),matches_neural_length=(b.ndim==1 and len(b)==len(t)))
   # Event counts from a published event table are descriptive only; new boundary
   # events need verified circle units/clock mapping and are not guessed here.
   if 'behavior/reversal_events' in f:
    event_shape=f['behavior/reversal_events'].shape
    reversal_count=int(event_shape[0]) if len(event_shape)==2 and event_shape[1]==2 else None
   else:reversal_count=None
   checks=['observed_signal_provenance_verified','behavior_clock_alignment_verified','timestamps_seconds_verified','identity_dictionary_join_verified','boundary_units_verified','session_boundaries_verified']
   unverified=[k for k in checks if not review.get(k,False)]
   rows.append(dict(animal=animal,neurons=len(x),frames=len(t),timestamps_strictly_increasing=time_ok,median_interval_native=float(np.median(dt)) if time_ok else None,gaps=int((dt>1.5*np.median(dt)).sum()) if time_ok else None,behavior=behavior,reversal_events_published_table_count=reversal_count,unverified_provenance_checks=unverified,usable_independent_animal_verified=False,reason='Numerical inspection does not by itself verify signal provenance, independent units or clock mapping'))
 save(LOCAL/'recording_audit_PRIVATE.json',rows);save(LOCAL/'neuron_quality_PRIVATE.json',nr);save(LOCAL/'recording_schemas_PRIVATE.json',schemas);save(LOCAL/'acquired_recording_hashes_PRIVATE.json',files)
 all_classes=sorted(set().union(*classes.values()) if classes else [])
 coverage=np.array([[c in classes[a] for a in [p.name.removesuffix('-cudata.h5') for p in paths]] for c in all_classes],bool)
 with (LOCAL/'neuron_class_by_animal_PRIVATE.csv').open('w') as out:
  w=csv.writer(out);w.writerow(['neuron_class']+[p.name.removesuffix('-cudata.h5') for p in paths])
  for c,v in zip(all_classes,coverage):w.writerow([c]+v.astype(int).tolist())
 with (LOCAL/'individual_neuron_by_animal_PRIVATE.csv').open('w') as out:
  w=csv.writer(out);animals=[p.name.removesuffix('-cudata.h5') for p in paths];w.writerow(['individual_neuron']+animals)
  for label in sorted(set().union(*individuals.values()) if individuals else []):w.writerow([label]+[int(label in individuals[a]) for a in animals])
 # Never declare readiness on finite values or supplied review flags alone.
 summary=dict(recordings_inspected=len(rows),genuine_usable_animal_count=None,neuron_classes_observed_in_metadata=len(all_classes),common_15class_15recording_numerical_gate=common_core_exists(coverage) if coverage.size else None,identity_dictionary_received=bool(labels),physiological_baselines_fitted=False,behavioral_event_alignment_verified=False,status='BLOCKED_NO_RECORDINGS' if not paths else 'SCIENTIFIC_PROVENANCE_REVIEW_REQUIRED',missing_evidence=['actual physiological payloads'] if not paths else ['manual confirmation of recording/identity provenance, units, independent animals and behavioural event alignment'])
 save(HERE/'RECORDING_QC_STATUS.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
