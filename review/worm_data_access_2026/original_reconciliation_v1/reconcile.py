"""Reconcile measured NWB target traces without fitting or interpolation.
Writes event-level observations only under ignored data/. Frozen QC unchanged.
"""
import argparse,csv,hashlib,importlib.util,json
from collections import Counter
from pathlib import Path
import sys
import numpy as np
import h5py

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=ROOT/'data/worm_data_access_2026/original_reconciliation_v1'
sys.path.insert(0,str(HERE.parent/'recording_qc_v1'))
import run_qc as frozen


def max_run(mask):
    mask=np.asarray(mask,bool)
    if not mask.any():return 0
    changes=np.diff(np.r_[False,mask,False].astype(int));return int(np.max(np.flatnonzero(changes==-1)-np.flatnonzero(changes==1)))


def finite_summary(v):
    v=np.asarray(v);good=np.isfinite(v)
    return dict(frames=len(v),finite=int(good.sum()),nan=int(np.isnan(v).sum()),infinite=int(np.isinf(v).sum()),
        finite_fraction=float(good.mean()) if len(v) else None,longest_nan_run_frames=max_run(~good))


def source_gate(t,y,start,lo,hi,fraction):
    # This is ONLY a run-length component on timestamp-selected NWB frames.
    # Original event volume indices/Dt and original nan_mask/matchless inputs are unavailable.
    sel=(t>=start+lo)&(t<start+hi);v=y[sel];ss=finite_summary(v)
    return dict(**ss,run_limit_frames=float(fraction*len(v)),nan_component_pass=bool(len(v)>0 and ss['longest_nan_run_frames']<=fraction*len(v)),
                clipped=bool(start+lo<t[0] or start+hi>t[-1]))


def main():
    p=argparse.ArgumentParser();p.add_argument('--wt-archive',type=Path);args=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);manifest=json.loads((HERE.parent/'access_resolution/selected_files.json').read_text())
    allrows=[];animals=[];schemas={}
    for entry in manifest:
        path=ROOT/'data/worm_nwb'/(entry['subject']+'.nwb');assert frozen.prepared.file_hash(path)==entry['sha256']
        with h5py.File(path,'r') as f:
            rt='processing/ophys/';g=f[rt+'GreenSignals/BaseGreenSignal'];red=f[rt+'RedSignals/BaseRedSignal'];interp=f[rt+'GreenSignals/InterpolatedGreenSignal']
            pp=f[rt+'PumpProbeGreenSegmentations/PumpProbeGreenPlaneSegmentation'];rp=f[rt+'PumpProbeRedSegmentations/PumpProbeRedPlaneSegmentation'];npal=f[rt+'NeuroPALSegmentations/NeuroPALPlaneSegmentation']
            names,conf,high,idcol=frozen.identities(pp,npal,g,.9)
            y=g['data'][:];t=g['timestamps'][:];rr=red['data'][:];it=interp['data'][:]
            st=f['intervals/OptogeneticStimulusTable'];starts=st['start_time'][:];stops=st['stop_time'][:];targets=st['target_pumpprobe_id'][:]
            assert np.array_equal(t,red['timestamps'][:]) and np.array_equal(t,interp['timestamps'][:])
            assert np.isfinite(t).all() and np.all(np.diff(t)>0) and np.max(np.diff(t))<3*np.median(np.diff(t))
            assert np.array_equal(pp['id'][:],rp['id'][:]) and np.array_equal(pp['neuropal_ids'][:],rp['neuropal_ids'][:])
            # Explicitly validate red as well as green column references; never assume ordering.
            region=red['table_region'];assert f[region.attrs['table']].name==rp.name
            assert np.array_equal(region[:],g['table_region'][:])
            qcs=frozen.trace_qc(t,y,starts,stops)
            # No field in these files stores original response/exclusion flags.
            paths=[]
            f.visititems(lambda n,o:paths.append(n) if isinstance(o,h5py.Dataset) else None)
            schemas[entry['subject']]=paths
            for i,v in enumerate(targets):
                if not np.isfinite(v) or v not in idcol or not high[idcol[v]]:continue
                col=idcol[v];start=starts[i];pre=(t>=start-10)&(t<start);post=(t>=start+1)&(t<=start+20)
                ps=finite_summary(y[pre,col]);qs=finite_summary(y[post,col]);target=y[:,col]
                obs_pre=y[pre,col][np.isfinite(y[pre,col])];obs_post=y[post,col][np.isfinite(y[post,col])]
                baseline=float(np.median(obs_pre)) if len(obs_pre) else None
                # Descriptive, observed-sample summary only. Not an autoresponse detector or statistical evidence.
                measurable=len(obs_pre)>=5 and len(obs_post)>=5 and baseline is not None and baseline>0
                delta=float((np.median(obs_post)-baseline)/baseline) if measurable else None
                af=source_gate(t,target,start,-30,60,.05);cf=source_gate(t,target,start,-30,60,.30);pf=source_gate(t,target,start,-30,30,.05)
                interwindow=it[(t>=start-10)&(t<=start+20),col]
                row=dict(subject=entry['subject'],session_start=frozen.txt(f['session_start_time'][()]),
                    stable_session_key=entry['filename'].split('_ses-')[1].split('_desc')[0].replace('-','_'),
                    event_index=i,nwb_total_events=len(starts),nwb_total_frames=len(t),nearest_nwb_stimulus_frame=int(np.argmin(np.abs(t-start))),target_pumpprobe_id=int(v),trace_column=int(col),target_name=names[col],confidence=float(conf[col]),
                    target_roi_present=True,green_target_trace_present=True,red_target_trace_present=True,
                    declared_identity_join_verified=True,unique_target_neuropal_link=True,
                    stimulus_onset_seconds=float(start),missing_timestamps=False,
                    pre_frames=ps['frames'],pre_finite=ps['finite'],pre_nan=ps['nan'],pre_finite_fraction=ps['finite_fraction'],pre_longest_nan_frames=ps['longest_nan_run_frames'],
                    post_frames=qs['frames'],post_finite=qs['finite'],post_nan=qs['nan'],post_finite_fraction=qs['finite_fraction'],post_longest_nan_frames=qs['longest_nan_run_frames'],
                    target_green_red_missing_mask_equal=bool(np.array_equal(np.isnan(target),np.isnan(rr[:,col]))),
                    stored_interpolated_window_finite_fraction=float(np.isfinite(interwindow).mean()),
                    alternate_representation='Stored InterpolatedGreenSignal, same ROI and timestamps; missing base values remain unmeasured',
                    same_frozen_windows_5pct_run_component_pass=bool(ps['longest_nan_run_frames']<=.05*ps['frames'] and qs['longest_nan_run_frames']<=.05*qs['frames']),
                    frozen_complete_target_trace_pass=qcs[i][col]['usable'],frozen_failure_reasons=';'.join(qcs[i][col]['reasons']),
                    frozen_baseline_f0=qcs[i][col]['f0'],frozen_baseline_noise=qcs[i][col]['noise'],
                    negative_samples_in_frozen_windows=qcs[i][col]['negative_samples'],zero_samples_in_frozen_windows=qcs[i][col]['zero_samples'],
                    nwb_60s_5pct_run_component_pass=pf['nan_component_pass'],nwb_90s_5pct_run_component_pass=af['nan_component_pass'],
                    nwb_90s_30pct_run_component_pass=cf['nan_component_pass'],
                    nwb_90s_longest_nan_run_frames=cf['longest_nan_run_frames'],nwb_90s_frames=cf['frames'],nwb_90s_window_clipped=cf['clipped'],
                    observed_target_has_ge5_pre_and_post_and_positive_baseline=measurable,
                    observed_only_median_post_relative_change=delta,
                    original_export_session_match='NOT_AVAILABLE',original_author_exclusion='UNKNOWN',
                    original_autoresponse_flag='UNKNOWN',original_published_inclusion='NOT_EVALUABLE',
                    missing_value_physical_cause='Unresolved: shared channel missingness is consistent with ROI extraction/tracking loss; original masks and logs needed',
                    artefact_cause='Not determined; no missing timestamp or negative-value cause in these target windows')
                allrows.append(row)
            animalrows=[r for r in allrows if r['subject']==entry['subject']]
            animals.append(dict(subject=entry['subject'],high_confidence_events=len(animalrows),
                target_rois_found=sum(r['target_roi_present'] for r in animalrows),
                any_missing_timestamps=sum(r['missing_timestamps'] for r in animalrows),
                failed_95pct_completeness=sum('finite_fraction_below_95pct' in r['frozen_failure_reasons'] for r in animalrows),
                failed_only_95pct_completeness=sum(r['frozen_failure_reasons']=='finite_fraction_below_95pct' for r in animalrows),
                same_frozen_windows_5pct_run_component_pass=sum(r['same_frozen_windows_5pct_run_component_pass'] for r in animalrows),
                long_missing_run_over_2s=sum('missing_run_over_2s' in r['frozen_failure_reasons'] for r in animalrows),
                low_or_undefined_snr=sum('low_baseline_snr' in r['frozen_failure_reasons'] for r in animalrows),
                nonpositive_or_absent_baseline=sum('nonpositive_baseline' in r['frozen_failure_reasons'] for r in animalrows),
                nwb_60s_5pct_run_component_pass=sum(r['nwb_60s_5pct_run_component_pass'] for r in animalrows),
                nwb_90s_5pct_run_component_pass=sum(r['nwb_90s_5pct_run_component_pass'] for r in animalrows),
                nwb_90s_30pct_run_component_pass=sum(r['nwb_90s_30pct_run_component_pass'] for r in animalrows),
                observed_target_summary_available=sum(r['observed_target_has_ge5_pre_and_post_and_positive_baseline'] for r in animalrows),
                full_trace_green_red_missing_masks_equal=bool(np.array_equal(np.isnan(y),np.isnan(rr))),
                observed_green_values_preserved_by_stored_interpolation=bool(np.array_equal(y[np.isfinite(y)],it[np.isfinite(y)])),
                original_full_inclusion_decisions_known=0))
    assert len(allrows)==29 and [r['high_confidence_events'] for r in animals]==[2,12,15]
    # Optional original text-archive reconciliation: stable session names, never OSF folder ordinal.
    original_status='NOT_AVAILABLE'
    if args.wt_archive:
        from inspect_exports import reconcile_text_archive
        original_status=reconcile_text_archive(args.wt_archive,allrows,manifest,OUT)
    with (OUT/'event_failure_decomposition.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=sorted({k for row in allrows for k in row}));w.writeheader();w.writerows(allrows)
    frozen.save(OUT/'nwb_dataset_inventory.json',schemas)
    with (HERE/'aggregate_failure_decomposition.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(animals[0]));w.writeheader();w.writerows(animals)
    # All row-level observations remain ignored; only across-event aggregates are published.
    sums={k:sum(r[k] for r in animals) for k in ['high_confidence_events','target_rois_found','any_missing_timestamps','failed_95pct_completeness','failed_only_95pct_completeness','same_frozen_windows_5pct_run_component_pass','long_missing_run_over_2s','low_or_undefined_snr','nonpositive_or_absent_baseline','nwb_60s_5pct_run_component_pass','nwb_90s_5pct_run_component_pass','nwb_90s_30pct_run_component_pass','observed_target_summary_available','original_full_inclusion_decisions_known']}
    summary=dict(nwb_reconciliation=sums,original_export_status=original_status,expansion_performed=False,
        source_component_checks_not_full_author_inclusion=True,frozen_qc_modified=False,no_interpolation_or_response_fit_performed=True,
        event_level_csv='data/worm_data_access_2026/original_reconciliation_v1/event_failure_decomposition.csv')
    frozen.save(HERE/'aggregate_summary.json',summary);print(json.dumps(summary,indent=2));print(json.dumps(animals,indent=2))
if __name__=='__main__':main()
