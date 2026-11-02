"""Frozen-rule full-cohort QC; no learning. Individual outputs exclusively private."""
from pathlib import Path
import json,hashlib,sys,tomllib,collections,csv
import numpy as np,pandas as pd,h5py
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];REVIEW=HERE.parent
BASE=ROOT/'data/worm_behavior_sheaf_readiness_v1';PRIVATE=BASE/'full_cohort_qc_v5'
sys.path[:0]=[str(REVIEW),str(REVIEW/'export_recovery_v3'),str(REVIEW/'sample_ingestion_v2')]
from qc_primitives import valid_windows,fixed_split
from reconcile_native_identities import primary_identities,common_core_witness
from sample_numerical_audit import exported_reversal_events,validate_transfer
from audit_behaviour_events import band_events
from extract_processed import digest

def clock_matches(tc,tn,sx,sy,xc,yc):
    """Same fixed tolerances and 1..16 contiguous-NIR bins as the sample audit."""
    offset=tn[0];ct=np.r_[0,np.cumsum(tn-offset)];cx=np.r_[0,np.cumsum(sx)];cy=np.r_[0,np.cumsum(sy)]
    mid=np.searchsorted(tn,tc);counts=np.zeros(len(tc),int)
    for shift in range(17):
        start=mid-shift
        for n in range(1,17):
            stop=start+n;ok=(start>=0)&(stop<=len(tn))&(start<=mid)&(mid<=stop)
            ix=np.flatnonzero(ok);a=start[ix];b=stop[ix]
            hit=(np.abs((ct[b]-ct[a])/n-(tc[ix]-offset))<=2e-7)&(np.abs((cx[b]-cx[a])/n-xc[ix])<2e-5)&(np.abs((cy[b]-cy[a])/n-yc[ix])<2e-5)
            counts[ix[hit]]+=1
    return {'frames':len(tc),'unique':int((counts==1).sum()),'ambiguous':int((counts>1).sum()),'unmatched':int((counts==0).sum())}

def plateau_lengths(x):
    """Exact repeated values are descriptive, not a sensor saturation diagnosis."""
    lengths=[]
    for row in x:
        breaks=np.r_[0,np.flatnonzero(np.diff(row)!=0)+1,len(row)]
        lengths.append(int(np.diff(breaks).max()))
    return lengths

def main():
    if (PRIVATE/'cohort_results_PRIVATE.json').exists():raise ValueError('Preserve existing full-cohort results; use versioned rerun')
    protocol=json.loads((REVIEW/'READINESS_PROTOCOL.json').read_text());frozen=json.loads((REVIEW/'FROZEN_PROTOCOL_SHA256.json').read_text())
    assert all(digest(ROOT/p)==h for p,h in frozen.items()),'Frozen file changed'
    extraction=json.loads((PRIVATE/'extraction_receipt_PRIVATE.json').read_text())
    for item in extraction['files']:
        p=BASE/'recordings'/item['filename'];assert p.stat().st_size==item['bytes'] and digest(p)==item['sha256'],'Extraction hash mismatch'
    assert validate_transfer()==5
    paths=sorted((BASE/'recordings').glob('*-cudata.h5'));animals=[p.name.removesuffix('-cudata.h5') for p in paths]
    expected={r['file'].removesuffix('-cudata.h5') for r in csv.DictReader((REVIEW/'PUBLIC_FILE_MANIFEST.csv').open()) if r['file'].endswith('-cudata.h5')}
    assert len(paths)==27 and set(animals)==expected
    splits=json.loads((BASE/'animal_splits_PRIVATE.json').read_text());assert splits==fixed_split(animals),'Split mismatch'
    native=BASE/'sample_ingestion_v2/native_identity';receipt=tomllib.loads((native/'native_receipt_PRIVATE.toml').read_text());table=native/'labels_PRIVATE.tsv'
    assert digest(table)==receipt['labels_sha256'] and digest(BASE/'recordings/NeuroPAL_labels.jld2')==receipt['source_sha256']
    labels=pd.read_csv(table,sep='\t',keep_default_na=False);assert set(labels.animal)==set(animals)
    primary=primary_identities(labels,protocol['identity_primary_confidence']);allrows=[];neuronrows=[];events=[];schemas={}
    for p,animal in zip(paths,animals):
        with h5py.File(p,'r') as f:
            schema={};f.visititems(lambda k,v:schema.update({k:{'shape':list(v.shape),'dtype':str(v.dtype)}}) if isinstance(v,h5py.Dataset) else None);schemas[animal]=schema
            t=f['timing/timestamp_confocal'][:];tn=f['timing/timestamp_nir'][:];x=f['gcamp/trace_array_original'][:].T;dt=np.diff(t)
            assert x.shape[1]==len(t) and np.isfinite(t).all() and (dt>0).all() and np.isfinite(tn).all() and (np.diff(tn)>0).all()
            segments=f['gcamp/idx_splits'][:].T
            assert np.array_equal(segments,[[1,len(t)]]),'Unknown segment layout; investigate instead of guessing'
            ids=labels[labels.animal==animal];hi=primary[primary.animal==animal]
            assert ids.trace_row_1based.between(1,len(x)).all() and not hi.trace_row_1based.duplicated().any() and not hi.label.duplicated().any()
            maps=[]
            for k in ['match_org_to_skip','match_skip_to_org']:
                m=f['gcamp/'+k][:].T;maps.append(bool(np.array_equal(m[:,0],m[:,1]) and set(m[:,0])==set(range(1,len(x)+1))))
            assert all(maps),'Unexpected ROI remapping requires explicit reconciliation'
            norms={k:f['gcamp/'+k][:].T for k in ['trace_array','traces_array_F_F20','traces_array_F_Fmean']}
            normalisation={'zscore_ddof1':bool(np.allclose(norms['trace_array'],(x-x.mean(1,keepdims=True))/x.std(1,ddof=1,keepdims=True),rtol=1e-10,atol=1e-10,equal_nan=True)), 'F_Fmean':bool(np.allclose(norms['traces_array_F_Fmean'],x/x.mean(1,keepdims=True),rtol=1e-10,atol=1e-10,equal_nan=True)), 'F_F20':bool(np.allclose(norms['traces_array_F_F20'],x/np.quantile(x,.2,axis=1,keepdims=True),rtol=1e-10,atol=1e-10,equal_nan=True))}
            b={k:f['behavior/'+k][:] for k in ['velocity','reversal_vec','head_angle','pumping','angular_velocity','dist_from_cuCenter']}
            assert all(len(v)==len(t) for v in b.values())
            rev=b['reversal_vec'];stored=f['behavior/reversal_events'][:].T
            assert np.array_equal(stored,exported_reversal_events(rev)),'Unexpected reversal convention'
            clocks=clock_matches(t,tn,f['behavior/stage_x_nir'][:],f['behavior/stage_y_nir'][:],f['behavior/stage_x'][:],f['behavior/stage_y'][:])
            plateaus=plateau_lengths(x);eligible_rows=[]
            for i,signal in enumerate(x):
                finite=float(np.isfinite(signal).mean());windows={str(h):int(valid_windows(signal,t,h).sum()) for h in [5,10,15]};constant=bool(np.nanstd(signal)==0)
                numerical=finite>=protocol['minimum_finite_trace_fraction'] and not constant and windows['10']>=protocol['minimum_valid_short_history_windows_per_eligible_neuron']
                match=hi[hi.trace_row_1based==i+1];identity=match.iloc[0] if len(match) else None
                row={'animal':animal,'trace_row_1based':i+1,'split':splits[animal],'finite_fraction':finite,'constant':constant,'numerical_windows':windows,'numerically_usable':bool(numerical),'primary_identity':identity is not None,'eligible_identified_processed_state':bool(numerical and identity is not None),'individual':str(identity.label) if identity is not None else '', 'class':str(identity.neuron_class) if identity is not None else '', 'LR':str(identity.LR) if identity is not None else '', 'DV':str(identity.DV) if identity is not None else '', 'negative_fraction':float((signal<0).mean()),'quantiles':np.nanquantile(signal,[0,.01,.25,.5,.75,.99,1]).tolist(),'max_exact_constant_run':plateaus[i],'genuine_observation_mask':None}
                neuronrows.append(row)
                if row['eligible_identified_processed_state']:eligible_rows.append(row)
            radius=float(f['behavior/circle_radius'][()])
            for band in [.1,.2,.3]:
                ev=band_events(t,b['dist_from_cuCenter'],radius,band,b['velocity'],rev)
                events.append({'animal':animal,'split':splits[animal],'band_mm':band,'events':ev})
            usable=bool(eligible_rows and all(np.isfinite(v).all() for v in b.values()) and clocks['unmatched']==0)
            allrows.append({'animal':animal,'split':splits[animal],'sha256':digest(p),'ROIs':len(x),'frames':len(t),'duration_seconds':float(t[-1]-t[0]),'dt_median':float(np.median(dt)),'dt_min':float(dt.min()),'dt_max':float(dt.max()),'neural_gaps':int((dt>1.5*np.median(dt)).sum()),'NIR_gaps':int((np.diff(tn)>1.5*np.median(np.diff(tn))).sum()),'all_original_finite':bool(np.isfinite(x).all()),'all_normalised_finite':all(np.isfinite(v).all() for v in norms.values()),'normalisation_relationships':normalisation,'clock_match':clocks,'primary_identified_individuals':len(hi),'eligible_individuals':len(eligible_rows),'eligible_classes':len({r['class'] for r in eligible_rows}),'side_resolved_eligible_individuals':sum(r['LR'] in ['L','R'] for r in eligible_rows),'processed_state_usable':usable,'reversal_intervals':len(stored),'initial_left_censored_reversal':bool(rev[0]),'subsequent_reversal_onsets':int((np.diff(rev)==1).sum()),'behaviour':{k:{'finite_fraction':float(np.isfinite(v).mean()),'std':float(np.nanstd(v))} for k,v in b.items()},'circle_radius_source_units':radius,'observation_provenance':'NO_ORIGINAL_MASK_OR_CANDIDATE_SPECIFIC_PIPELINE_CONFIG'})
    private={'animals':allrows,'neurons':neuronrows,'events':events,'schemas':schemas,'split':splits}
    (PRIVATE/'cohort_results_PRIVATE.json').write_text(json.dumps(private,indent=2)+'\n')
    nd=pd.DataFrame([{k:v for k,v in row.items() if k not in ['quantiles','numerical_windows']} for row in neuronrows]);eligible=nd[nd.eligible_identified_processed_state]
    usable=[r['animal'] for r in allrows if r['processed_state_usable']]
    coverage={}
    for key,name in [('class','classes'),('individual','individuals')]:
        c=pd.crosstab(eligible[key],eligible.animal).reindex(columns=animals,fill_value=0)>0
        c.astype(int).to_csv(PRIVATE/f'{name}_by_animal_PRIVATE.csv');w=common_core_witness(c[usable].values)
        coverage[name]={'distinct':len(c),'common_all27':int(c.all(axis=1).sum()),'same15_in_same15_usable_animals':w is not None,'identities_in_at_least15_animals':int((c.sum(axis=1)>=15).sum())}
        if w is not None:(PRIVATE/f'{name}_common_core_witness_PRIVATE.json').write_text(json.dumps({'animals':[usable[i] for i in w[0]],'identities':[c.index[i] for i in w[1]],'QC_only_not_feature_selection':True},indent=2)+'\n')
    circuit=[]
    for label in protocol['diagnostic_circuit_individual_neurons']:
        a=set(eligible[eligible.individual==label].animal)
        circuit.append({'individual_neuron':label,'all_eligible_animals':len(a),**{s+'_animals':sum(splits[i]==s for i in a) for s in ['train','validation','test']}})
    pd.DataFrame(circuit).to_csv(HERE/'FROZEN_CIRCUIT_COVERAGE_AGGREGATE.csv',index=False)
    ea=[]
    for band in [.1,.2,.3]:
        selected=[r for r in events if r['band_mm']==band];ev=[e for r in selected for e in r['events']];risk=[e for e in ev if e['at_risk_for_new_reversal']]
        ea.append({'band_mm':band,'animals':len(selected),'animals_with_events':sum(bool(r['events']) for r in selected),'animals_with_repeated_events':sum(len(r['events'])>=2 for r in selected),'total_events':len(ev),'event_count_range':[min(len(r['events']) for r in selected),max(len(r['events']) for r in selected)],'secondary_at_risk':len(risk),'secondary_positives':sum(e['new_reversal_within3s'] for e in risk),'velocity_target_std':float(np.std([e['velocity_target'] for e in ev])) if ev else None,'by_split':{s:{'animals_with_events':sum(bool(r['events']) for r in selected if r['split']==s),'events':sum(len(r['events']) for r in selected if r['split']==s),'secondary_at_risk':sum(e['at_risk_for_new_reversal'] for r in selected if r['split']==s for e in r['events']),'secondary_positives':sum(e['new_reversal_within3s'] and e['at_risk_for_new_reversal'] for r in selected if r['split']==s for e in r['events'])} for s in ['train','validation','test']},'independently_calibrated_head_geometry':False,'stimulus_TTL':False})
    for summary in ea:
        grouped=[r for r in events if r['band_mm']==summary['band_mm']]
        summary['animals_with_at_risk_event']=sum(any(x['at_risk_for_new_reversal'] for x in row['events']) for row in grouped)
        summary['animals_with_secondary_positive']=sum(any(x['at_risk_for_new_reversal'] and x['new_reversal_within3s'] for x in row['events']) for row in grouped)
        summary['secondary_negatives']=summary['secondary_at_risk']-summary['secondary_positives']
    (HERE/'BEHAVIOURAL_QC_AGGREGATE.json').write_text(json.dumps(ea,indent=2)+'\n')
    agg={'actual_independent_animals':len(paths),'single_sessions':len(paths),'technical_processed_state_usable_animals':len(usable),'usable_animals_by_frozen_split':dict(collections.Counter(splits[a] for a in usable)),'total_ROI_recordings':len(neuronrows),'ROI_count_range':[min(r['ROIs'] for r in allrows),max(r['ROIs'] for r in allrows)],'frames_per_recording':sorted({r['frames'] for r in allrows}),'duration_range_seconds':[min(r['duration_seconds'] for r in allrows),max(r['duration_seconds'] for r in allrows)],'median_sampling_interval_range_seconds':[min(r['dt_median'] for r in allrows),max(r['dt_median'] for r in allrows)],'minimum_maximum_sampling_interval_seconds':[min(r['dt_min'] for r in allrows),max(r['dt_max'] for r in allrows)],'all_original_and_normalised_finite':all(r['all_original_finite'] and r['all_normalised_finite'] for r in allrows),'all_normalisation_relationships_verified':all(all(r['normalisation_relationships'].values()) for r in allrows),'constant_ROI_recordings':sum(r['constant'] for r in neuronrows),'neural_gap_count':sum(r['neural_gaps'] for r in allrows),'NIR_gap_count':sum(r['NIR_gaps'] for r in allrows),'numerically_usable_ROI_recordings':sum(r['numerically_usable'] for r in neuronrows),'eligible_identified_individual_recordings':len(eligible),'eligible_individuals_per_animal_range':[min(r['eligible_individuals'] for r in allrows),max(r['eligible_individuals'] for r in allrows)],'eligible_classes_per_animal_range':[min(r['eligible_classes'] for r in allrows),max(r['eligible_classes'] for r in allrows)],'side_resolved_eligible_recordings':sum(r['side_resolved_eligible_individuals'] for r in allrows),'coverage':coverage,'valid_numerical_windows_by_history':{str(h):sum(r['numerical_windows'][str(h)] for r in neuronrows) for h in [5,10,15]},'clock_frames':sum(r['clock_match']['frames'] for r in allrows),'clock_unique_frames':sum(r['clock_match']['unique'] for r in allrows),'clock_ambiguous_frames':sum(r['clock_match']['ambiguous'] for r in allrows),'clock_unmatched_frames':sum(r['clock_match']['unmatched'] for r in allrows),'reversal_intervals':sum(r['reversal_intervals'] for r in allrows),'initial_left_censored_animals':sum(r['initial_left_censored_reversal'] for r in allrows),'subsequent_reversal_onsets':sum(r['subsequent_reversal_onsets'] for r in allrows),'max_exact_constant_run_samples':max(r['max_exact_constant_run'] for r in neuronrows),'ROIs_with_negative_samples':sum(r['negative_fraction']>0 for r in neuronrows),'original_observation_masks_available':False,'saturation_detector_metadata_available':False,'candidate_specific_preprocessing_available':False,'frozen_protocol_hashes_unchanged':True,'model_fitted':False}
    dates=collections.defaultdict(set)
    for animal,split in splits.items():dates[animal[:10]].add(split)
    agg.update(distinct_recording_SHA256=len({r['sha256'] for r in allrows}),distinct_acquisition_dates=len(dates),dates_spanning_frozen_splits=sum(len(groups)>1 for groups in dates.values()),exported_dataset_counts=sorted({len(schema) for schema in schemas.values()}))
    (HERE/'COHORT_QC_AGGREGATE.json').write_text(json.dumps(agg,indent=2)+'\n');print(json.dumps(agg,indent=2));print(json.dumps(ea,indent=2))
if __name__=='__main__':main()
