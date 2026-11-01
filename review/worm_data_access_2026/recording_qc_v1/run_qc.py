"""Offline recording QC only; details/arrays written exclusively to ignored data/.
Rules frozen in QC_RULES.json. No response, dynamical or sheaf model fitting.
"""
import csv
import hashlib
import importlib.util
import io
import json
from collections import Counter
from itertools import combinations
from pathlib import Path
import sys
import h5py
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
LOCAL=ROOT/'data/worm_data_access_2026/recording_qc_v1'
DATA=ROOT/'data/worm_nwb'
RULES=json.loads((HERE/'QC_RULES.json').read_text())
spec=importlib.util.spec_from_file_location('prepared',HERE.parent/'access_resolution/audit_transferred.py')
prepared=importlib.util.module_from_spec(spec);spec.loader.exec_module(prepared)


def txt(x):return x.decode() if isinstance(x,bytes) else str(x)
def save(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def envtext(name):return json.loads((ROOT/'data/worm_data_access_2026'/f'{name}.json').read_text())['structuredContent']['content']
CANON={l.split()[1] for l in envtext('neuron_ids').splitlines()}


def longest_missing(t,finite,dt):
    idx=np.flatnonzero(~finite)
    if not len(idx):return 0.
    runs=np.split(idx,np.flatnonzero(np.diff(idx)>1)+1)
    return float(max(t[r[-1]]-t[r[0]]+dt for r in runs))


def trace_qc(t,y,starts,stops,min_snr=5):
    """Return per-event x ROI masks/reasons. No amplitude-response selection."""
    dt=float(np.median(np.diff(t)));out=[];flat=np.zeros(y.shape[1],bool)
    for j in range(y.shape[1]):
        z=y[np.isfinite(y[:,j]),j];flat[j]=len(z)<2 or np.ptp(z)==0
    for i,start in enumerate(starts):
        pre=(t>=start-10)&(t<start);post=(t>=start+1)&(t<=start+20)
        flags=[]
        if start-10<t[0] or start+20>t[-1]:flags.append('window_boundary')
        if np.any((np.arange(len(starts))!=i)&(starts<start+20)&(stops>start-10)):flags.append('other_stimulus_overlap')
        selected=(t>=start-10)&(t<=start+20)
        if np.any(np.diff(t[selected])>3*dt):flags.append('timestamp_gap')
        rr=[]
        for j in range(y.shape[1]):
            why=list(flags);p=y[pre,j];q=y[post,j];pf=np.isfinite(p);qf=np.isfinite(q)
            if min(len(p),len(q))<5:why.append('too_few_frames')
            if not len(p) or not len(q) or pf.mean()<.95 or qf.mean()<.95:why.append('finite_fraction_below_95pct')
            if longest_missing(t[pre],pf,dt)>2 or longest_missing(t[post],qf,dt)>2:why.append('missing_run_over_2s')
            f0=float(np.median(p[pf])) if pf.any() else None
            adj=pf[:-1]&pf[1:];dif=np.diff(p)[adj]
            noise=float(np.median(np.abs(dif-np.median(dif)))/.67448975/np.sqrt(2)) if len(dif) else None
            if f0 is None or f0<=0:why.append('nonpositive_baseline')
            elif noise is None or f0<min_snr*noise:why.append('low_baseline_snr')
            if flat[j]:why.append('flat_or_unobserved_recording_roi')
            v=np.r_[p[pf],q[qf]]
            rr.append(dict(usable=not why,reasons=sorted(set(why)),f0=f0,noise=noise,
                negative_samples=int((v<0).sum()),zero_samples=int((v==0).sum())))
        out.append(rr)
    return out


def identities(pp,npal,signal,min_conf):
    reg=signal['table_region'];rows=reg[:].astype(int)
    assert signal.file[reg.attrs['table']].name==pp.name
    assert len(rows)==signal['data'].shape[1] and len(set(rows))==len(rows)
    pids=pp['id'][:].tolist();nids=npal['id'][:].tolist();labels=[txt(x).strip() for x in npal['labels'][:]]
    conf=npal['labels_confidences'][:];mapping={k:(lab,float(c)) for k,lab,c in zip(nids,labels,conf)}
    links=[txt(x).strip() for x in pp['neuropal_ids'][:]]
    named=[];confidence=[]
    for row in rows:
        link=links[row];lab,c=mapping.get(int(link),('',-1)) if link.isdigit() else ('',-1)
        named.append(lab);confidence.append(c)
    ct=Counter(x for x in named if x)
    eligible=np.array([lab in CANON and ct[lab]==1 and np.isfinite(c) and min_conf<=c<=1 for lab,c in zip(named,confidence)])
    return named,np.asarray(confidence),eligible,{pids[row]:col for col,row in enumerate(rows)}


def audit_one(e):
    path=DATA/(e['subject']+'.nwb');assert path.stat().st_size==e['bytes']
    sha=prepared.file_hash(path);assert sha==e['sha256']
    with h5py.File(path,'r') as f:
        rt='processing/ophys/';g=f[rt+'GreenSignals/BaseGreenSignal'];y=g['data'][:];t=g['timestamps'][:]
        st=f['intervals/OptogeneticStimulusTable'];starts=st['start_time'][:];stops=st['stop_time'][:];target=st['target_pumpprobe_id'][:]
        pp=f[rt+'PumpProbeGreenSegmentations/PumpProbeGreenPlaneSegmentation'];npal=f[rt+'NeuroPALSegmentations/NeuroPALPlaneSegmentation']
        names,conf,high,id_to_col=identities(pp,npal,g,RULES['primary_identity_minimum'])
        good=trace_qc(t,y,starts,stops);good_snr={str(v):trace_qc(t,y,starts,stops,v) for v in RULES['baseline_snr_sensitivity']}
        metrics=prepared.array_metrics(t,y,starts,stops)
        manual=[];event_details=[];receiver_masks=[];edges=[]
        priorities=['missing_manual_target','target_has_no_recorded_trace','target_not_confident_unique_canonical','target_trace_fails_qc','no_usable_other_high_confidence_receiver']
        reasons=Counter();validtargets=0;accepted=0;targetqcp=0
        for i,val in enumerate(target):
            why=[];col=None
            if not np.isfinite(val) or not float(val).is_integer() or val<0:why.append(priorities[0])
            elif val not in id_to_col:why.append(priorities[1])
            else:
                col=id_to_col[val]
                if not high[col]:why.append(priorities[2])
                else:validtargets+=1
            if col is not None and not good[i][col]['usable']:why.append(priorities[3])
            if col is not None and high[col] and good[i][col]['usable']:targetqcp+=1
            rec=np.array([high[j] and j!=col and good[i][j]['usable'] for j in range(len(names))])
            if not rec.any():why.append(priorities[4])
            ok=not why;accepted+=int(ok);receiver_masks.append(rec)
            if why:reasons[why[0]]+=1
            if ok:
                for j in np.flatnonzero(rec):edges.append((names[col],names[j]))
            event_details.append(dict(event_index=i,start_s=float(starts[i]),stop_s=float(stops[i]),
                stored_manual_target_id=float(val) if np.isfinite(val) else None,target_name=names[col] if col is not None else None,
                target_confidence=float(conf[col]) if col is not None else None,accepted=ok,
                reasons=why,usable_receivers=int(rec.sum()),target_trace_qc=good[i][col] if col is not None else None))
        finite=y[np.isfinite(y)];neg=(y<0)&np.isfinite(y);negative_cols=np.flatnonzero(neg.any(axis=0))
        neg_windows=sum(any(x['negative_samples'] for x in rr) for rr in good)
        target_identity_sensitivity={}
        for threshold in [1.,.9,.8]:
            _,_,hi,_=identities(pp,npal,g,threshold)
            target_identity_sensitivity[str(threshold)]=sum(bool(np.isfinite(v) and v in id_to_col and hi[id_to_col[v]]) for v in target)
        num_all_pairs=0;receiver_reason=Counter();high_pairtotal=0
        for i,rr in enumerate(good):
            for j,q in enumerate(rr):
                if high[j]:
                    high_pairtotal+=1;num_all_pairs+=q['usable']
                    for reason in q['reasons']:receiver_reason[reason]+=1
        sensitivity={}
        for key,mask in [('no_negative_in_target_or_receiver',None),('no_zero_in_target_or_receiver',None),*[(f'baseline_snr_{v}',z) for v,z in good_snr.items()]]:
            n=0
            for i,ev in enumerate(event_details):
                if ev['target_name'] is None:continue
                val=target[i];col=id_to_col.get(val)
                if col is None or not high[col]:continue
                rr=good[i] if mask is None else mask[i]
                field='negative_samples' if key.startswith('no_negative') else 'zero_samples' if key.startswith('no_zero') else None
                if not rr[col]['usable'] or (field and rr[col][field]):continue
                if any(high[j] and j!=col and rr[j]['usable'] and (not field or not rr[j][field]) for j in range(len(names))):n+=1
            sensitivity[key]=n
        repeated=list(set(tuple(x) for x in edges))
        named_set={x for x in names if x in CANON}
        base_interp=f[rt+'GreenSignals/InterpolatedGreenSignal'];red=f[rt+'RedSignals/BaseRedSignal']
        sid=txt(f['general/subject/subject_id'][()]);assert sid==e['subject'][4:]
        rawtargetvalid=sum(np.isfinite(v) and v in id_to_col for v in target)
        r=dict(subject=e['subject'],hash_verified=True,sha256=sha,session_start=txt(f['session_start_time'][()]),
            filename=path.name,**metrics,manually_assigned_recorded_target_events=int(rawtargetvalid),
            canonical_named_rois=sum(x in CANON for x in names),unique_canonical_names=len(named_set),
            high_confidence_unique_rois=int(high.sum()),high_confidence_manual_target_events=validtargets,
            high_confidence_manual_target_fraction=validtargets/len(starts),target_qc_passing_events=targetqcp,
            accepted_events=accepted,rejected_events=len(starts)-accepted,acceptance_fraction=accepted/len(starts),
            exclusive_rejection_reasons=dict(reasons),identity_threshold_sensitivity=target_identity_sensitivity,
            event_qc_sensitivity=sensitivity,events_with_ge5_receivers=sum(ev['accepted'] and ev['usable_receivers']>=5 for ev in event_details),
            events_with_ge10_receivers=sum(ev['accepted'] and ev['usable_receivers']>=10 for ev in event_details),
            high_confidence_receiver_event_pairs=high_pairtotal,qc_passing_high_confidence_receiver_event_pairs=int(num_all_pairs),
            overlapping_receiver_trace_rejection_reasons=dict(receiver_reason),
            fluorescence_quantiles=dict(zip(['min','p001','p01','median','p99','max'],np.quantile(finite,[0,.001,.01,.5,.99,1]).tolist())),
            negative_finite_samples=int(neg.sum()),negative_fraction_of_finite=float(neg.sum()/len(finite)),
            negative_rois=len(negative_cols),negative_roi_names=[names[j] for j in negative_cols],
            zero_fraction_of_finite=float((finite==0).mean()),trial_windows_touching_any_negative_roi=neg_windows,
            negative_within_accepted_receiver_trials=sum(good[i][j]['negative_samples']>0 for i,ev in enumerate(event_details) if ev['accepted'] for j in np.flatnonzero(receiver_masks[i])),
            cross_channel_timestamps_equal=bool(np.array_equal(t,red['timestamps'][:]) and np.array_equal(t,base_interp['timestamps'][:])),
            interpolated_nan_fraction=float(np.isnan(base_interp['data'][:]).mean()),
            start_time_s=float(t[0]),stop_time_s=float(t[-1]),stimulus_start_range_s=[float(starts.min()),float(starts.max())],
            minimum_interstimulus_interval_s=float(np.diff(starts).min()),median_interstimulus_interval_s=float(np.median(np.diff(starts))),
            stored_power_values=sorted(set(st['power'][:].tolist())),unit=txt(g['data'].attrs['unit']),
            target_field_description=txt(st['target_pumpprobe_id'].attrs['description']),
            high_confidence_names=sorted(set(np.array(names)[high])),canonical_named_names=sorted(named_set),
            high_confidence_stimulated_names=sorted({ev['target_name'] for ev in event_details if ev['target_name'] and ev['target_confidence'] is not None and ev['target_confidence']>=.9 and ev['target_name'] in set(np.array(names)[high])}),
            accepted_stimulated_names=sorted({ev['target_name'] for ev in event_details if ev['accepted']}),
            accepted_sender_receiver_pairs=sorted(repeated),events=event_details,
            per_roi=[dict(column=j,name=names[j],confidence=float(conf[j]),high_confidence=bool(high[j]),
                 finite_fraction=float(np.isfinite(y[:,j]).mean()),negative_samples=int(neg[:,j].sum()),
                 min=float(np.min(y[np.isfinite(y[:,j]),j])) if np.isfinite(y[:,j]).any() else None) for j in range(y.shape[1])])
        # Keep detailed receiver masks and reason tables local, not in report directory.
        save(LOCAL/(e['subject']+'_details.json'),r)
        save(LOCAL/(e['subject']+'_receiver_qc.json'),good)
        return r


def anatomy_audit(records):
    rows=list(csv.DictReader(io.StringIO(envtext('connectome')),delimiter='\t'))
    chemical={(r['pre'],r['post']) for r in rows if r['type']=='chemical' and r['pre'] in CANON and r['post'] in CANON}
    gap={tuple(sorted([r['pre'],r['post']])) for r in rows if r['type']=='electrical' and r['pre'] in CANON and r['post'] in CANON}
    allsets=[set(tuple(e) for e in r['accepted_sender_receiver_pairs']) for r in records]
    counts=Counter(e for s in allsets for e in s);results=[]
    for edge,count in sorted(counts.items(),key=lambda x:(-x[1],x[0])):
        results.append(dict(sender=edge[0],receiver=edge[1],animals_with_eligible_trial=count,
            chemical_direct_in_witvliet7=edge in chemical,electrical_direct_in_witvliet7=tuple(sorted(edge)) in gap))
    save(LOCAL/'anatomical_pair_support.json',results)
    return dict(reference='Witvliet specimen 7; reference animal not recorded animals; chemical rows directed, electrical canonicalised undirected',
        candidate_sender_receiver_pairs=len(results),candidate_pairs_in_ge2_animals=sum(r['animals_with_eligible_trial']>=2 for r in results),
        direct_chemical_candidates=sum(r['chemical_direct_in_witvliet7'] for r in results),
        direct_electrical_candidates=sum(r['electrical_direct_in_witvliet7'] for r in results),
        repeated_direct_candidates=[r for r in results if r['animals_with_eligible_trial']>=2 and (r['chemical_direct_in_witvliet7'] or r['electrical_direct_in_witvliet7'])],
        repeated_effective_candidates_without_direct_reference_edge=[r for r in results if r['animals_with_eligible_trial']>=2 and not (r['chemical_direct_in_witvliet7'] or r['electrical_direct_in_witvliet7'])])


def main():
    LOCAL.mkdir(parents=True,exist_ok=True)
    entries=json.loads((HERE.parent/'access_resolution/selected_files.json').read_text())
    records=[audit_one(e) for e in entries]
    overlap={}
    for a,b in combinations(records,2):
        overlap[a['subject']+'__'+b['subject']]=dict(shared_high_confidence_recorded_names=len(set(a['high_confidence_names'])&set(b['high_confidence_names'])),
            shared_high_confidence_manual_targets=sorted(set(a['high_confidence_stimulated_names'])&set(b['high_confidence_stimulated_names'])),
            shared_accepted_targets=sorted(set(a['accepted_stimulated_names'])&set(b['accepted_stimulated_names'])))
    inter=set.intersection(*(set(r['high_confidence_names']) for r in records))
    anatomy=anatomy_audit(records)
    summary=dict(animals=len(records),total_events=sum(r['trials'] for r in records),
        high_confidence_manual_target_events=sum(r['high_confidence_manual_target_events'] for r in records),
        accepted_events=sum(r['accepted_events'] for r in records),
        shared_high_confidence_recorded_names_all_animals=sorted(inter),pairwise=overlap,anatomy=anatomy,
        schema_correction='Actual NWBs use table_region, not rois. HDF5 region reference and order explicitly verified.',
        versions=dict(python=sys.version.split()[0],numpy=np.__version__,h5py=h5py.__version__,hdf5=h5py.version.hdf5_version),
        rules_sha256=hashlib.sha256((HERE/'QC_RULES.json').read_bytes()).hexdigest())
    save(LOCAL/'summary.json',summary)
    # Requested per-animal aggregate QC only; no traces, event identities or masks in review/.
    fields=['subject','frames','rois','canonical_named_rois','unique_canonical_names','high_confidence_unique_rois','trials','manually_assigned_recorded_target_events','high_confidence_manual_target_events','accepted_events','rejected_events','acceptance_fraction','nan_fraction','negative_finite_samples','negative_rois']
    with (HERE/'animal_qc_counts.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:r[k] for k in fields} for r in records)
    print(json.dumps({k:v for k,v in summary.items() if k!='anatomy'},indent=2))
    print('Anatomy summary',json.dumps(anatomy,indent=2))
    for r in records:print('Animal QC',json.dumps({k:r[k] for k in fields},indent=2),r['exclusive_rejection_reasons'])
if __name__=='__main__':main()
