"""Actual small-sample HDF5 audit. Layout correction, not a gate amendment.
No identity heuristic, fitted neural model, imputation or reconstruction of traces.
All animal/neuron results remain in ignored data; public results aggregate only.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np,h5py
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
PRIVATE=ROOT/'data/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2'
PARENT=HERE.parent;sys.path.insert(0,str(PARENT));from qc_primitives import valid_windows

def canonical_traces(matrix,n_times):
 m=np.asarray(matrix)
 if m.ndim!=2 or m.shape[0]!=n_times:raise ValueError('Not the verified Julia HDF5 time-by-neuron storage layout')
 return m.T

def canonical_reversal_table(matrix):
 m=np.asarray(matrix)
 if m.ndim!=2 or m.shape[0]!=2:raise ValueError('Not the verified Julia 2-by-events storage layout')
 return m.T

def exported_reversal_events(binary):
 r=np.asarray(binary)
 if not np.isin(r,[0,1]).all():raise ValueError('Nonbinary reversal state')
 starts=np.flatnonzero(np.diff(r)==1)+1;ends=np.flatnonzero(np.diff(r)==-1)+1
 if r[0]:starts=np.r_[1,starts] # Reproduce author's mixed initial/start convention, not a new onset estimate.
 if len(starts)!=len(ends):raise ValueError('Censored final reversal requires separate interpretation')
 return np.c_[starts,ends]

def validate_transfer():
 base=PRIVATE.parent;receipt=json.loads((base/'transfer_receipt_PRIVATE.json').read_text(encoding='utf-8-sig'))
 if receipt['algorithm'].upper().replace('-','')!='SHA256' or len(receipt['files'])!=5:raise ValueError('Unexpected receipt')
 names=set()
 for row in receipt['files']:
  name=row['filename']
  if Path(name).name!=name or name in names:raise ValueError('Invalid/duplicate receipt filename')
  names.add(name);p=base/'recordings'/name
  if p.stat().st_size!=int(row['bytes']) or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256'].lower():raise ValueError('Transfer mismatch')
 return len(names)

def main():
 for name,digest in json.loads((PARENT/'FROZEN_PROTOCOL_SHA256.json').read_text()).items():
  if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Frozen protocol changed')
 count=validate_transfer();rows=[];neurons=[]
 for p in sorted((PRIVATE.parent/'recordings').glob('*-cudata.h5')):
  with h5py.File(p,'r') as f:
   t=f['timing/timestamp_confocal'][:];tn=f['timing/timestamp_nir'][:];dt=np.diff(t);dn=np.diff(tn)
   if not np.isfinite(t).all() or np.any(dt<=0):raise ValueError('Invalid neural timestamps')
   x=canonical_traces(f['gcamp/trace_array_original'][:],len(t))
   norms={key:canonical_traces(f['gcamp/'+key][:],len(t)) for key in ['trace_array','traces_array_F_F20','traces_array_F_Fmean']}
   expected_z=(x-x.mean(1,keepdims=True))/x.std(1,ddof=1,keepdims=True)
   expected_mean=x/x.mean(1,keepdims=True);expected_q=x/np.quantile(x,.2,axis=1,keepdims=True)
   segment_pairs=f['gcamp/idx_splits'][:].T
   if not np.array_equal(segment_pairs,[[1,len(t)]]):raise ValueError('New session layout requires source review')
   for i,trace in enumerate(x):
    neurons.append(dict(animal=p.stem,trace_row_1based=i+1,finite_fraction=float(np.isfinite(trace).mean()),constant=bool(np.nanstd(trace)==0),valid_numerical_history_windows=int(valid_windows(trace,t,10).sum()),negative_fraction=float((trace<0).mean()),minimum=float(np.nanmin(trace)),maximum=float(np.nanmax(trace)),observed_versus_imputed_unknown=True))
   bstats={}
   for key in ['velocity','reversal_vec','head_angle','pumping','angular_velocity','dist_from_cuCenter']:
    b=f['behavior/'+key][:]
    bstats[key]=dict(n=len(b),finite_fraction=float(np.isfinite(b).mean()),std=float(np.nanstd(b)))
   rev=f['behavior/reversal_vec'][:];events=canonical_reversal_table(f['behavior/reversal_events'][:])
   event_match=np.array_equal(events,exported_reversal_events(rev))
   if not event_match:raise ValueError('Stored reversal events differ from the inspected author exporter')
   maps=[]
   for key in ['match_org_to_skip','match_skip_to_org']:
    m=f['gcamp/'+key][:].T;maps.append(bool(np.array_equal(m[:,0],m[:,1]) and set(m[:,0])==set(range(1,len(x)+1))))
   # Descriptive only: do not certify circles/head geometry or genuine observations.
   row=dict(animal=p.name.removesuffix('-cudata.h5'),neurons=len(x),frames=len(t),duration_seconds=float(t[-1]-t[0]),median_neural_interval_seconds=float(np.median(dt)),min_neural_interval_seconds=float(dt.min()),max_neural_interval_seconds=float(dt.max()),neural_gaps_over_1_5_median=int((dt>1.5*np.median(dt)).sum()),nir_frames=len(tn),median_nir_interval_seconds=float(np.median(dn)),nir_gaps_over_1_5_median=int((dn>1.5*np.median(dn)).sum()),all_raw_finite=bool(np.isfinite(x).all()),all_norm_finite=all(np.isfinite(v).all() for v in norms.values()),zscore_matches_ddof1=bool(np.allclose(norms['trace_array'],expected_z,rtol=1e-10,atol=1e-10)),mean_normalisation_matches=bool(np.allclose(norms['traces_array_F_Fmean'],expected_mean,rtol=1e-10,atol=1e-10)),q20_normalisation_matches=bool(np.allclose(norms['traces_array_F_F20'],expected_q,rtol=1e-10,atol=1e-10)),single_recording_segment=True,row_mapping_is_identity_relation=all(maps),stored_reversal_events=len(events),reversal_table_reproduced=True,initial_reversal_left_censored=bool(rev[0]),behavior=bstats,neuron_identity_QC='PENDING_NATIVE_JLD2',original_observation_masks='NOT PRESENT IN 30 EXPORTED DATASETS',circle_geometry_reconciliation='NOT VERIFIED: stage-origin restoration alone does not reproduce supplied radial distances')
   rows.append(row)
 (PRIVATE/'sample_numeric_results_PRIVATE.json').write_text(json.dumps(rows,indent=2)+'\n');(PRIVATE/'neuron_numeric_results_PRIVATE.json').write_text(json.dumps(neurons,indent=2)+'\n')
 clocks=json.loads((PRIVATE/'reconstructed_clock_bins_PRIVATE.json').read_text())
 agg=dict(transfer_integrity='PASS',windows_vm_hashes_and_sizes_matching=count,recordings_numerically_inspected=len(rows),neuron_count_range=[min(r['neurons'] for r in rows),max(r['neurons'] for r in rows)],total_neuron_recordings=sum(r['neurons'] for r in rows),independent_animals_not_neurons=len(rows),frames_per_recording=sorted(set(r['frames'] for r in rows)),all_raw_and_normalised_values_finite=all(r['all_raw_finite'] and r['all_norm_finite'] for r in rows),constant_neuron_recordings=sum(n['constant'] for n in neurons),neural_timestamp_gaps_over_1_5_median=sum(r['neural_gaps_over_1_5_median'] for r in rows),valid_numerical_10history_windows=sum(n['valid_numerical_history_windows'] for n in neurons),normalisation_relationships_verified=all(r['zscore_matches_ddof1'] and r['mean_normalisation_matches'] and r['q20_normalisation_matches'] for r in rows),reversal_event_count_range=[min(r['stored_reversal_events'] for r in rows),max(r['stored_reversal_events'] for r in rows)],total_reversal_intervals=sum(r['stored_reversal_events'] for r in rows),initial_left_censored_reversal_animals=sum(r['initial_reversal_left_censored'] for r in rows),clock_frames_jointly_matched=sum(r['summary']['unique_joint_bins']+r['summary']['ambiguous_bins'] for r in clocks),clock_frames_with_unique_joint_match=sum(r['summary']['unique_joint_bins'] for r in clocks),clock_ambiguous_frames=sum(r['summary']['ambiguous_bins'] for r in clocks),identity_status='NOT EVALUABLE until native Julia/JLD2',genuine_observation_fraction=None,full_sample_technical_gate='PENDING_IDENTITY_AND_OBSERVATION_PROVENANCE',full_27_acquisition_started=False,model_fitted=False)
 (HERE/'SAMPLE_QC_AGGREGATE.json').write_text(json.dumps(agg,indent=2)+'\n');print(json.dumps(agg,indent=2))
if __name__=='__main__':main()
