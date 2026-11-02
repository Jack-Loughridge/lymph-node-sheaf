"""Aggregate-only cluster uncertainty and complete pair/model reporting."""
from pathlib import Path
import json,hashlib,collections
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];PRIVATE=ROOT/'data/worm_behavior_forecasting_v1'
CONTRASTS=(('VAR','AR'),('VAR_beh','AR_beh'),('rank1','AR'),('graph','AR'),('AR_beh','AR'),('AR','persistence'),('VAR','persistence'),('VAR_beh','persistence'))
B=10000

def stable_seed(token):return (20261102+int(hashlib.sha256(token.encode()).hexdigest()[:8],16))%(2**32)

def contrast_stats(base,full,groups=None,seed=20261102):
    base=np.asarray(base,float);full=np.asarray(full,float);n=len(base)
    if n<2 or not (base>0).all():raise ValueError('Insufficient or degenerate independent error pairs')
    groups=np.arange(n).astype(str) if groups is None else np.asarray(groups)
    unique=sorted(set(groups));a=np.array([base[groups==g].sum() for g in unique]);b=np.array([full[groups==g].sum() for g in unique])
    rng=np.random.default_rng(seed);draw=rng.integers(0,len(unique),(B,len(unique)))
    bootstrap=100*(1-b[draw].sum(1)/a[draw].sum(1))
    point=100*(1-full.mean()/base.mean());lo,hi=np.quantile(bootstrap,[.025,.975])
    jack=[100*(1-np.delete(full,i).mean()/np.delete(base,i).mean()) for i in range(n)]
    return {'relative_MSE_reduction_percent':float(point),'CI95_low_percent':float(lo),'CI95_high_percent':float(hi),'n_animals':n,'n_bootstrap_clusters':len(unique),'animals_improved':int((full<base).sum()),'animals_worsened':int((full>base).sum()),'median_animal_reduction_percent':float(np.median(100*(1-full/base))),'evaluation_leave_one_animal_out_min_percent':float(min(jack)),'evaluation_leave_one_animal_out_max_percent':float(max(jack)),'baseline_mean_MSE':float(base.mean()),'augmented_mean_MSE':float(full.mean())}

def main():
    results=json.loads((PRIVATE/'all_results_PRIVATE.json').read_text());flat=[];modelrows=[];contrasts=[];status=[]
    for result in results:
        common={'neuron_1':result['pair'][0],'neuron_2':result['pair'][1],'history':result['history'],'grouping':result['grouping']}
        status.append({**common,'status':result['status'],'coobserved_animals':result['independent_animals'],'distinct_grouping_units':result['distinct_groups'],'evaluated_animal_folds':len(result['folds'])})
        if result['status']!='OK':continue
        for receiver in range(2):
            identity={'sender':result['pair'][1-receiver],'receiver':result['pair'][receiver],'history':result['history'],'grouping':result['grouping']}
            for fold in result['folds']:
                for model,metrics in fold['metrics'].items():
                    flat.append({**identity,'animal':fold['animal'],'day':fold['acquisition_day'],'model':model,'MSE_raw':metrics['MSE_raw'][receiver],'MSE_normalised':metrics['MSE_normalised'][receiver],'training_MSE_raw':metrics['training_MSE_raw'][receiver],'training_MSE_normalised':metrics['training_MSE_normalised'][receiver],'alpha':metrics['alpha'][receiver],'converged':metrics['converged'][receiver],'n_windows':fold['n_windows']})
    frame=pd.DataFrame(flat);frame.to_csv(PRIVATE/'fold_metrics_PRIVATE.csv',index=False)
    for key,g in frame.groupby(['sender','receiver','history','grouping'],sort=True):
        ident=dict(zip(['sender','receiver','history','grouping'],key))
        for model,m in g.groupby('model',sort=True):
            modelrows.append({**ident,'model':model,'n_animals':len(m),'n_days':m.day.nunique(),'total_target_windows':int(m.n_windows.sum()),'outer_mean_MSE_raw':m.MSE_raw.mean(),'outer_mean_MSE_normalised':m.MSE_normalised.mean(),'apparent_training_mean_MSE_raw':m.training_MSE_raw.mean(),'apparent_training_mean_MSE_normalised':m.training_MSE_normalised.mean(),'alpha_outer_fold_counts':json.dumps(dict(sorted(collections.Counter(m.alpha.astype(str)).items()))),'outer_nonconverged_predictions':int((~m.converged).sum())})
        for full,base in CONTRASTS:
            a=g[g.model==base].set_index('animal').sort_index();b=g[g.model==full].set_index('animal').reindex(a.index)
            assert len(a)==len(b) and a.n_windows.equals(b.n_windows)
            for scale in ['normalised','raw']:
                stats=contrast_stats(a['MSE_'+scale],b['MSE_'+scale],seed=stable_seed(str(key)+scale))
                row={**ident,'comparison':full+'_vs_'+base,'scale':scale,'interval_type':'animal_percentile_conditional_on_CV_predictions',**stats}
                if key[-1]=='day':
                    clustered=contrast_stats(a['MSE_'+scale].to_numpy(),b['MSE_'+scale].to_numpy(),a.day.to_numpy(),stable_seed(str(key)+scale+'date'))
                    row.update(day_cluster_CI95_low_percent=clustered['CI95_low_percent'],day_cluster_CI95_high_percent=clustered['CI95_high_percent'],day_clusters=clustered['n_bootstrap_clusters'])
                contrasts.append(row)
    pd.DataFrame(modelrows).to_csv(HERE/'MODEL_COMPARISONS.csv',index=False);cf=pd.DataFrame(contrasts);cf.to_csv(HERE/'DIRECTED_CONTRASTS.csv',index=False)
    pd.DataFrame(status).to_csv(HERE/'PAIR_EVALUABILITY.csv',index=False)
    # Average directions within an animal first; never count overlapping pair targets as animals.
    pooled=[]
    for (history,grouping),g in frame.groupby(['history','grouping']):
        for full,base in CONTRASTS:
            a=g[g.model==base].groupby('animal').MSE_normalised.mean().sort_index();b=g[g.model==full].groupby('animal').MSE_normalised.mean().reindex(a.index)
            stats=contrast_stats(a,b,seed=stable_seed(str((history,grouping,full,base))))
            row={'history':history,'grouping':grouping,'comparison':full+'_vs_'+base,'scale':'normalised','aggregation':'mean supported directions within animal then equal animal weights',**stats}
            if grouping=='day':
                day=np.array([x[:10] for x in a.index]);dc=contrast_stats(a.to_numpy(),b.to_numpy(),day,stable_seed(str((history,grouping,full,base))+'day'))
                row.update(day_cluster_CI95_low_percent=dc['CI95_low_percent'],day_cluster_CI95_high_percent=dc['CI95_high_percent'],day_clusters=dc['n_bootstrap_clusters'])
            pooled.append(row)
    pd.DataFrame(pooled).to_csv(HERE/'POOLED_CONTRASTS.csv',index=False)
    matched=[]
    for history in sorted(frame.history.unique()):
        day=frame[(frame.history==history)&(frame.grouping=='day')]
        keys=set(zip(day.sender,day.receiver,day.animal))
        for grouping in ['animal','day']:
            g=frame[(frame.history==history)&(frame.grouping==grouping)]
            g=g[[k in keys for k in zip(g.sender,g.receiver,g.animal)]]
            for full,base in CONTRASTS:
                aa=g[g.model==base].groupby('animal').MSE_normalised.mean().sort_index();bb=g[g.model==full].groupby('animal').MSE_normalised.mean().reindex(aa.index)
                stat=contrast_stats(aa,bb,seed=stable_seed(str((history,grouping,full,base,'matched'))))
                matched.append({'history':history,'grouping':grouping,'comparison':full+'_vs_'+base,'scope':'same animal-direction support as evaluable day folds',**stat})
    pd.DataFrame(matched).to_csv(HERE/'DAY_MATCHED_POOL.csv',index=False)
    primary=cf[(cf.history==10)&(cf.grouping=='animal')&(cf.scale=='normalised')]
    summary={'eligible_unordered_pairs':10,'directed_relationships':20,'evaluated_training_animals_union':int(frame.animal.nunique()),'union_acquisition_days':int(frame.day.nunique()),'primary':{contrast:{'directions':len(s),'directions_positive_point_estimate':int((s.relative_MSE_reduction_percent>0).sum()),'directions_CI_excluding_zero_positive':int((s.CI95_low_percent>0).sum()),'directions_meeting_5percent_and_positive_CI':int(((s.relative_MSE_reduction_percent>=5)&(s.CI95_low_percent>0)).sum()),'directions_improve_majority_animals':int((s.animals_improved>s.n_animals/2).sum()),'median_direction_reduction_percent':float(s.relative_MSE_reduction_percent.median())} for contrast,s in primary.groupby('comparison')},'bootstrap_draws':B,'interval_caveats':['descriptive unadjusted direction intervals','3-6 coobserved animals per pair','conditional on CV predictions without refitting every bootstrap','overlapping CV training folds'],'validation_test_recordings_used':False,'full_sheaf_fitted':False}
    (HERE/'RESULTS_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2));print(pd.DataFrame(pooled).query("history==10 and comparison in ['VAR_vs_AR','VAR_beh_vs_AR_beh']").to_string(index=False))
if __name__=='__main__':main()
