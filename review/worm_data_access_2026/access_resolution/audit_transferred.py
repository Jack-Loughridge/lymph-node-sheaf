"""Offline audit of exactly the approved processed NWBs; no models or raw imaging.
Only stdlib is needed for --verify-only. Full inspection needs numpy + h5py.
Detailed reports contain animal-level metadata: save them only in ignored data/.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent


def file_hash(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def text(v):
    return v.decode('utf-8',errors='replace') if isinstance(v,bytes) else str(v)


def array_metrics(t,y,starts,stops):
    import numpy as np
    t=np.asarray(t,float);y=np.asarray(y,float);starts=np.asarray(starts,float);stops=np.asarray(stops,float)
    if t.ndim!=1 or len(t)<2 or y.ndim!=2 or y.shape[0]!=len(t):raise ValueError('Expected time x ROI signals with matching timestamps')
    if not np.isfinite(t).all() or np.any(np.diff(t)<=0):raise ValueError('Timestamps nonfinite or not strictly increasing')
    if starts.ndim!=1 or stops.shape!=starts.shape or not np.isfinite(starts).all() or not np.isfinite(stops).all() or np.any(stops<starts):raise ValueError('Malformed stimulation intervals')
    dt=np.diff(t);med=float(np.median(dt));finite=np.isfinite(y);counts=finite.sum(axis=0)
    nonconstant=0
    for j in range(y.shape[1]):
        v=y[finite[:,j],j]
        if len(v)>1 and np.ptp(v)>0:nonconstant+=1
    idx=np.clip(np.searchsorted(t,starts),0,len(t)-1);prev=np.maximum(idx-1,0)
    nearest=np.minimum(abs(t[idx]-starts),abs(t[prev]-starts))
    inside=(starts>=t[0]) & (stops<=t[-1])
    gaps=dt>3*med
    baseline=[];post=[]
    for start in starts:
        # A previous event's 30s response interval must not overlap the pre-window.
        clean_pre=start-30>=t[0] and not np.any((starts<start)&(stops+30>start-30))
        clean_pre=clean_pre and not np.any(gaps&(t[:-1]<start)&(t[1:]>start-30))
        clean_post=start+30<=t[-1] and not np.any((starts>start)&(starts<start+30))
        clean_post=clean_post and not np.any(gaps&(t[:-1]<start+30)&(t[1:]>start))
        baseline.append(bool(clean_pre));post.append(bool(clean_post))
    return dict(frames=len(t),rois=y.shape[1],trials=len(starts),duration_s=float(t[-1]-t[0]),
        median_dt_s=med,max_gap_s=float(dt.max()),gaps_gt_3_median=int(gaps.sum()),
        nan_fraction=float(np.isnan(y).mean()),infinite_values=int(np.isinf(y).sum()),
        rois_all_missing=int((counts==0).sum()),rois_with_nonconstant_finite_trace=nonconstant,
        roi_finite_fraction_quantiles=np.quantile(finite.mean(axis=0),[0,.25,.5,.75,1]).tolist(),
        events_outside_time_range=int((~inside).sum()),events_out_of_order=int((np.diff(starts)<0).sum()),
        event_start_nearest_frame_distance_s_quantiles=np.quantile(nearest,[0,.5,1]).tolist() if len(nearest) else [],
        events_more_than_one_median_dt_from_frame=int((nearest>med).sum()),
        stimulus_duration_s_values=sorted(set((stops-starts).tolist())),
        clean_30s_pre_windows_after_30s_previous_response_washout=int(sum(baseline)),
        gap_free_uncontaminated_30s_post_windows=int(sum(post)))


def inspect_nwb(path,entry):
    import numpy as np
    import h5py
    with h5py.File(path,'r') as f:
        root='processing/ophys/'
        g=f[root+'GreenSignals/BaseGreenSignal'];st=f['intervals/OptogeneticStimulusTable']
        pp=f[root+'PumpProbeGreenSegmentations/PumpProbeGreenPlaneSegmentation']
        npal=f[root+'NeuroPALSegmentations/NeuroPALPlaneSegmentation']
        t=g['timestamps'][:];y=g['data'][:];starts=st['start_time'][:];stops=st['stop_time'][:]
        r=array_metrics(t,y,starts,stops)
        warnings=[];errors=[]
        if r['infinite_values']:errors.append('Infinite signal values')
        if r['events_outside_time_range']:errors.append('Events outside calcium recording time range')
        if r['events_out_of_order']:errors.append('Stimuli not chronological')
        ids=npal['id'][:].tolist();labels=[text(x).strip() for x in npal['labels'][:]]
        if len(ids)!=len(set(ids)):raise ValueError('Nonunique NeuroPAL ids')
        if len(ids)!=len(labels):raise ValueError('NeuroPAL id-label lengths differ')
        lmap=dict(zip(ids,labels));confmap={}
        if 'labels_confidences' in npal:
            conf=npal['labels_confidences'][:]
            if len(conf)!=len(ids):raise ValueError('NeuroPAL confidence length mismatch')
            confmap=dict(zip(ids,[text(x) for x in conf]))
            r['confidence_dataset_description']=text(npal['labels_confidences'].attrs.get('description',''))
        else:warnings.append('No NeuroPAL confidence dataset')
        links=[text(x).strip() for x in pp['neuropal_ids'][:]];roi_ids=pp['id'][:].tolist()
        if len(links)!=len(roi_ids):raise ValueError('ROI id-link lengths differ')
        if len(set(roi_ids))!=len(roi_ids):raise ValueError('Nonunique PumpProbe ROI ids')
        if 'rois' in g:
            roi_rows=np.asarray(g['rois'][:],int)
            if 'table' in g['rois'].attrs:
                ref=g['rois'].attrs['table']
                if isinstance(ref,h5py.Reference) and ref and f[ref].name!=pp.name:raise ValueError('GreenSignal ROI table differs from expected PumpProbe table')
            else:warnings.append('GreenSignal ROI table reference missing; inspect column semantics')
        elif len(links)==y.shape[1]:
            roi_rows=np.arange(len(links));warnings.append('No signal ROIs table: positional identity mapping provisional')
        else:raise ValueError('Cannot associate trace columns with ROIs')
        if roi_rows.ndim!=1 or len(roi_rows)!=y.shape[1] or np.any(roi_rows<0) or np.any(roi_rows>=len(links)) or len(set(roi_rows.tolist()))!=len(roi_rows):raise ValueError('Invalid trace-column ROI table rows')
        column_labels=[];confidence_counts=Counter();noninteger_links=0;unmatched_numeric=0
        for row in roi_rows:
            link=links[row]
            if link.isdigit():
                nid=int(link);lab=lmap.get(nid,'')
                if nid not in lmap:unmatched_numeric+=1
                confidence_counts[confmap.get(nid,'MISSING')]+=1
            else:
                lab='';noninteger_links+=int(bool(link))
            column_labels.append(lab)
        usable=[x for x in column_labels if x and x.lower() not in {'none','nan','unknown','unidentified','?'}]
        r.update(labelled_roi_count=len(usable),labelled_roi_fraction=len(usable)/len(column_labels) if column_labels else None,
                 distinct_recorded_neuron_names=len(set(usable)),duplicate_nonempty_neuron_names=len(usable)-len(set(usable)),
                 labels_with_question_mark=sum('?' in x for x in column_labels),unmatched_numeric_neuropal_links=unmatched_numeric,
                 other_nonempty_neuropal_links=noninteger_links,confidence_value_histogram=dict(confidence_counts),
                 identity_reliability='Label presence measured; confidence conventions and independent corroboration still require review.')
        # The field names a PumpProbe id. Resolve ids, then signal-table rows; never nearest-soma matching.
        target_values=st['target_pumpprobe_id'][:]
        if target_values.shape!=starts.shape:raise ValueError('Target/event length mismatch')
        id_to_row={v:i for i,v in enumerate(roi_ids)};row_to_col={int(row):i for i,row in enumerate(roi_rows)}
        targeted_names=[];assigned=0;recorded=0;target_counts=Counter()
        for val in target_values:
            key=val.item() if hasattr(val,'item') else val
            target_counts[str(key)]+=1
            if isinstance(key,(int,float)) and np.isfinite(key) and float(key).is_integer() and key>=0 and key in id_to_row:
                assigned+=1;row=id_to_row[key]
                if row in row_to_col:
                    recorded+=1;lab=column_labels[row_to_col[row]]
                    if lab in usable:targeted_names.append(lab)
        r.update(assigned_target_trials=assigned,unresolved_target_trials=len(starts)-assigned,
                 target_recorded_in_green_trials=recorded,neuron_named_target_trials=len(targeted_names),
                 distinct_named_stimulated_neurons=len(set(targeted_names)),stored_target_value_histogram=dict(target_counts),
                 target_field_description=text(st['target_pumpprobe_id'].attrs.get('description','')),
                 target_mapping='target_pumpprobe_id -> PumpProbe table id -> GreenSignal rois row. No inferred targets.')
        if not np.array_equal(pp['id'][:],np.arange(len(links))):warnings.append('Nonsequential ROI ids: confirm target id versus positional-index semantics with conversion code')
        for group in ['GreenSignals/InterpolatedGreenSignal','RedSignals/BaseRedSignal']:
            key=group.replace('/','_')
            if root+group in f:
                alt=f[root+group];r[key+'_timestamps_equal']=bool(np.array_equal(t,alt['timestamps'][:]))
                r[key+'_shape']=list(alt['data'].shape)
                if not r[key+'_timestamps_equal']:warnings.append(group+' timestamps differ from base green')
            else:warnings.append(group+' absent')
        r['green_signal_unit']=text(g['data'].attrs.get('unit',''))
        if 'power' in st:
            power=st['power'][:];r['stored_power_values']=sorted(set(np.asarray(power).tolist()))
            r['stored_power_unit']=text(st['power'].attrs.get('unit',''))
            if len(power)!=len(starts):errors.append('Power column length mismatch')
        r['dose_caveat']='Stored duration/power can be conversion constants, not measured delivered dose.'
        sid=text(f['general/subject/subject_id'][()]);r['subject_id']=sid
        if sid.removeprefix('sub-')!=entry['subject'].removeprefix('sub-'):errors.append('Subject differs from selected manifest')
        for field in ['strain','genotype','age','growth_stage']:
            if 'general/subject/'+field in f:r[field]=text(f['general/subject/'+field][()])
        r['session_start']=text(f['session_start_time'][()])
        if 'centroids' in pp:
            c=pp['centroids'][:];r['centroid_shape']=list(c.shape);r['centroids_all_finite']=bool(np.isfinite(c).all())
            r['centroid_unit_attr']=text(pp['centroids'].attrs.get('unit','UNSPECIFIED'))
        if r['duplicate_nonempty_neuron_names']:warnings.append('Neuron names duplicated across recorded ROI columns')
        if r['gaps_gt_3_median']:warnings.append('Long timestamp gaps: do not interpolate silently')
        r.update(warnings=warnings,structural_errors=errors,structure_status='FAIL' if errors else 'STRUCTURE_CHECKS_PASS',
                 biological_validation='NOT PERFORMED',recorded_names=sorted(set(usable)),stimulated_names=sorted(set(targeted_names)))
        return r


def main():
    p=argparse.ArgumentParser();p.add_argument('--data-dir',type=Path,required=True);p.add_argument('--manifest',type=Path,default=HERE/'selected_files.json')
    p.add_argument('--out-dir',type=Path);p.add_argument('--verify-only',action='store_true');args=p.parse_args()
    entries=json.loads(args.manifest.read_text());out=args.out_dir or args.data_dir.parent/'transferred_sample_audit'
    if out.resolve()==HERE.resolve() or 'review' in out.resolve().parts:raise SystemExit('Save animal-level audit under ignored data/, not review/')
    out.mkdir(parents=True,exist_ok=True);results=[];fail=False
    for e in entries:
        path=args.data_dir/e['filename'];r={'filename':e['filename']}
        try:
            if '_desc-segmentation_' not in path.name:raise ValueError('Not a processed segmentation file')
            if not path.is_file():raise ValueError('Missing file')
            if path.stat().st_size!=e['bytes']:raise ValueError('Size mismatch')
            got=file_hash(path)
            if got!=e['sha256']:raise ValueError('SHA256 mismatch')
            with path.open('rb') as f:
                if f.read(8)!=b'\x89HDF\r\n\x1a\n':raise ValueError('HDF5 magic mismatch')
            r.update(integrity='PASS',sha256=got)
            if not args.verify_only:
                r.update(inspect_nwb(path,e))
                if r['structural_errors']:fail=True
        except Exception as ex:
            r.update(error=type(ex).__name__+': '+str(ex));fail=True
        results.append(r)
    (out/'recordings.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n')
    parsed=[r for r in results if 'subject_id' in r]
    counts=Counter(n for r in parsed for n in set(r['recorded_names']));stim=Counter(n for r in parsed for n in set(r['stimulated_names']))
    # This sample has one expected recording per animal. Detect duplicate animal/session/payload before counting.
    subjects={r['subject_id'] for r in parsed};sessions={(r['subject_id'],r['session_start']) for r in parsed}
    if len(subjects)!=len(parsed) and parsed:fail=True
    hashes=[r['sha256'] for r in parsed]
    if len(set(hashes))!=len(hashes):fail=True
    aggregate=dict(status='INCOMPLETE_OR_FAILED' if fail else ('CHECKSUMS_ONLY_PASS' if args.verify_only else 'STRUCTURE_CHECKS_PASS_REVIEW_QC'),
        expected_animals=len({e['subject'] for e in entries}),files_with_verified_checksums=sum(r.get('integrity')=='PASS' for r in results),
        recording_payloads_inspected=len(parsed),independent_subject_ids=len(subjects),distinct_subject_sessions=len(sessions),
        stimulation_trials=sum(r['trials'] for r in parsed) if parsed else None,
        recorded_neuron_names_observed_in_at_least_two_animals=sum(v>=2 for v in counts.values()) if parsed else None,
        stimulated_neuron_names_observed_in_at_least_two_animals=sum(v>=2 for v in stim.values()) if parsed else None,
        recorded_name_animal_coverage_histogram=dict(Counter(counts.values())),
        note='Integrity and naming audit only. Identity credibility, biological independence and predictability require review; no sheaf fit.')
    (out/'aggregate.json').write_text(json.dumps(aggregate,indent=2)+'\n')
    print(json.dumps(aggregate,indent=2));print('Detailed animal-level results saved locally under '+str(out))
    return 1 if fail else 0
if __name__=='__main__':sys.exit(main())
