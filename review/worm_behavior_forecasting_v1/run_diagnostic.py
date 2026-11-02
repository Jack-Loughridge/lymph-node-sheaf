"""Training-only nested animal/date forecasts. No sheaf, no held-out cohort access."""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import sys,json,hashlib,itertools,tomllib,datetime,csv
import numpy as np,pandas as pd,h5py
from models import MODELS,GRID,moments,combine,fit,mse_from_moments,predict
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
OLD=ROOT/'review/worm_behavior_sheaf_readiness_v1';BASE=ROOT/'data/worm_behavior_sheaf_readiness_v1';PRIVATE=ROOT/'data/worm_behavior_forecasting_v1'
sys.path[:0]=[str(OLD),str(OLD/'export_recovery_v3')]
from qc_primitives import valid_windows
from reconcile_native_identities import primary_identities

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def require_training(animal,split):
    if split.get(animal)!='train':raise ValueError('Validation/test or unknown animal access refused')

def designs(neural,beh,t,h):
    n=len(t);valid=np.ones(n-h,bool)
    for channel in np.r_[neural,beh]:valid &= valid_windows(channel,t,h)
    target=np.arange(h,n)[valid]
    # Explicit indices make the cutoff auditable: most recent feature is target−1.
    idx=target[:,None]-np.arange(1,h+1)[None,:]
    x=np.concatenate([channel[idx] for channel in np.r_[neural,beh]],axis=1)
    return x,neural[:,target].T,target

def source_data(pair,split,identities):
    eligible=[];data={};hashes={}
    for animal in sorted(a for a,s in split.items() if s=='train'):
        ids=identities[identities.animal==animal]
        if not all((ids.label==n).sum()==1 for n in pair):continue
        require_training(animal,split)
        path=BASE/'recordings'/f'{animal}-cudata.h5'
        with h5py.File(path,'r') as f:
            t=f['timing/timestamp_confocal'][:];rows=[int(ids[ids.label==n].trace_row_1based.iloc[0])-1 for n in pair]
            neural=np.stack([f['gcamp/trace_array_original'][:,i] for i in rows])
            assert np.array_equal(f['gcamp/idx_splits'][:].T,[[1,len(t)]])
            v=f['behavior/velocity'][:]
            beh=np.stack([v,np.abs(v),f['behavior/reversal_vec'][:],f['behavior/dist_from_cuCenter'][:],f['behavior/pumping'][:]])
            if not all(np.isfinite(c).mean()>=.95 and np.std(c)>0 and valid_windows(c,t,10).sum()>=200 for c in neural):continue
            data[animal]=(neural,beh,t);eligible.append(animal);hashes[animal]=sha(path)
    return eligible,data,hashes

def graph_support(pair,reference):
    result=[];details=[]
    for receiver in range(2):
        sender=1-receiver;a=pair[sender];b=pair[receiver]
        chem=reference[(reference.pre==a)&(reference.post==b)&(reference.type=='chemical')]
        elec=reference[(((reference.pre==a)&(reference.post==b))|((reference.pre==b)&(reference.post==a)))&(reference.type=='electrical')]
        result.append(bool(len(chem) or len(elec)))
        details.append({'sender':a,'receiver':b,'chemical_listed':bool(len(chem)),'electrical_listed_either_orientation':bool(len(elec)),'cross_history_allowed':result[-1]})
    return tuple(result),details

def nested(pair,h,kind,arrays,graph,runhash):
    destination=PRIVATE/'checkpoints'/f'{pair[0]}_{pair[1]}_H{h}_{kind}_PRIVATE.json'
    if destination.exists():
        result=json.loads(destination.read_text())
        if result['run_hash']!=runhash:raise ValueError('Checkpoint differs from frozen run; preserve and version')
        return result
    ids=sorted(arrays);groups={a:(a if kind=='animal' else a[:10]) for a in ids};ug=sorted(set(groups.values()))
    result={'pair':pair,'history':h,'grouping':kind,'run_hash':runhash,'independent_animals':len(ids),'distinct_groups':len(ug),'folds':[],'status':'OK','fit_count':0,'nonconverged_rank1_fits':0,'max_rank1_objective_increase':0.}
    if len(ug)<3:
        result['status']='NOT EVALUABLE: fewer than three independent grouping units'
        destination.write_text(json.dumps(result,indent=2)+'\n');return result
    ms={a:moments(arrays[a][0],arrays[a][1]) for a in ids};cache={};scale_cache={}
    def cached(train,alpha,model):
        key=(tuple(sorted(train)),float(alpha),model)
        if key not in cache:
            skey=key[0]
            if skey not in scale_cache:scale_cache[skey]=combine([ms[a] for a in skey])
            f=scale_cache[skey];b,offset,info=fit(f,h,alpha,model,graph)
            cache[key]=(f,b,offset,info);result['fit_count']+=1
            if model=='rank1':
                result['nonconverged_rank1_fits']+=not info['converged']
                result['max_rank1_objective_increase']=max(result['max_rank1_objective_increase'],info['maximum_objective_increase'])
        return cache[key]
    for outer in ug:
        held=[a for a in ids if groups[a]==outer];train=[a for a in ids if groups[a]!=outer];innergroups=sorted({groups[a] for a in train})
        selected={};innerlog={}
        for model in MODELS:
            if model=='persistence':selected[model]=[0.,0.];continue
            scores=[]
            for alpha in GRID:
                losses=[]
                for inner in innergroups:
                    fitids=[a for a in train if groups[a]!=inner];valid=[a for a in train if groups[a]==inner]
                    f,b,o,info=cached(fitids,alpha,model)
                    losses.extend(mse_from_moments(ms[a],f,b,o)[1].tolist() for a in valid)
                scores.append(np.mean(losses,axis=0))
            scores=np.array(scores);innerlog[model]={'alphas':GRID,'mean_animal_normalised_MSE':scores.tolist()}
            if model=='rank1':
                alpha=min(range(len(GRID)),key=lambda i:(scores[i].mean(),-GRID[i]));selected[model]=[GRID[alpha]]*2
            else:selected[model]=[GRID[min(range(len(GRID)),key=lambda i:(scores[i,j],-GRID[i]))] for j in range(2)]
        # Both output directions use exactly the same windows and outer information cutoff.
        for animal in held:
            predictions={};metrics={};y=arrays[animal][1];target=arrays[animal][2]
            for model in MODELS:
                pred=np.empty_like(y);raw=[];normal=[];traw=[];tnorm=[];conv=[]
                for j in range(2):
                    f,b,o,info=cached(train,selected[model][j],model)
                    pred[:,j]=predict(arrays[animal][0],f,b,o)[:,j]
                    loss=float(np.mean((pred[:,j]-y[:,j])**2));raw.append(loss);normal.append(loss/f.sy[j]**2)
                    trainloss=[mse_from_moments(ms[a],f,b,o) for a in train]
                    traw.append(float(np.mean([r[0][j] for r in trainloss])));tnorm.append(float(np.mean([r[1][j] for r in trainloss])))
                    conv.append(info['converged'])
                predictions[model]=pred;metrics[model]={'MSE_raw':raw,'MSE_normalised':normal,'training_MSE_raw':traw,'training_MSE_normalised':tnorm,'alpha':selected[model],'converged':conv}
            predfile=PRIVATE/'predictions'/f'{pair[0]}_{pair[1]}_H{h}_{kind}_{animal}_PRIVATE.npz'
            if predfile.exists():raise ValueError('Uncheckpointed prediction already exists; preserve and inspect')
            np.savez_compressed(predfile,targets=y,target_indices=target,**predictions)
            result['folds'].append({'animal':animal,'acquisition_day':animal[:10],'outer_group':outer,'fit_animals':train,'n_windows':len(y),'metrics':metrics,'inner_scores':innerlog,'prediction_sha256':sha(predfile)})
        # Checkpoint after a complete panel, avoiding any mixing with modified code.
    destination.write_text(json.dumps(result,indent=2)+'\n');return result

def main():
    PRIVATE.mkdir(exist_ok=True)
    for sub in ['checkpoints','predictions']:(PRIVATE/sub).mkdir(exist_ok=True)
    prereg=json.loads((HERE/'PREREGISTRATION_RECEIPT.json').read_text())
    assert all(sha(ROOT/p)==h for p,h in prereg['reference_SHA256'].items()),'Preregistered source changed'
    frozen=json.loads((OLD/'FROZEN_PROTOCOL_SHA256.json').read_text());assert all(sha(ROOT/p)==h for p,h in frozen.items())
    config=json.loads((OLD/'READINESS_PROTOCOL.json').read_text());split=json.loads((BASE/'animal_splits_PRIVATE.json').read_text());assert sorted(list(split.values()).count(s) for s in ['train','validation','test'])==[6,6,15]
    native=BASE/'sample_ingestion_v2/native_identity';nr=tomllib.loads((native/'native_receipt_PRIVATE.toml').read_text());assert sha(native/'labels_PRIVATE.tsv')==nr['labels_sha256']
    labels=primary_identities(pd.read_csv(native/'labels_PRIVATE.tsv',sep='\t',keep_default_na=False));labels=labels[labels.animal.isin([a for a,s in split.items() if s=='train'])]
    reference=pd.read_csv(ROOT/'data/bioelectric_focused_feasibility_2026/white_whole.tsv',sep='\t')
    sources={str(p.relative_to(ROOT)):sha(p) for p in [HERE/'models.py',HERE/'run_diagnostic.py',HERE/'IMPLEMENTATION_PREREGISTRATION.md',OLD/'READINESS_PROTOCOL.json']}
    runhash=hashlib.sha256(json.dumps(sources,sort_keys=True).encode()).hexdigest()
    if (PRIVATE/'run_receipt_PRIVATE.json').exists():
        assert json.loads((PRIVATE/'run_receipt_PRIVATE.json').read_text())['run_hash']==runhash,'Use versioned run for changed code'
    receipt={'run_hash':runhash,'code_and_protocol':sources,'started_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'allowed_split':'train','recording_sources_PRIVATE':{},'validation_test_recordings_opened':False}
    (PRIVATE/'run_receipt_PRIVATE.json').write_text(json.dumps(receipt,indent=2)+'\n')
    expected=pd.read_csv(OLD/'full_cohort_qc_v5/TRAINING_PAIR_PANEL_COVERAGE.csv');results=[];graphrows=[]
    for pair in itertools.combinations(config['diagnostic_circuit_individual_neurons'],2):
        ids,data,hashes=source_data(pair,split,labels)
        # Guard eligibility against a silent change to the frozen support table.
        old=expected[(expected.neuron_1==pair[0])&(expected.neuron_2==pair[1])].iloc[0]
        assert len(ids)==int(old.training_coobserved_animals)
        receipt['recording_sources_PRIVATE'].update(hashes)
        (PRIVATE/'run_receipt_PRIVATE.json').write_text(json.dumps(receipt,indent=2)+'\n')
        if len(ids)<3:continue
        graph,details=graph_support(pair,reference);graphrows.extend(details)
        for h in [10,5,15]:
            arrays={a:designs(*data[a],h) for a in ids}
            assert all(len(values[1])>=200 for values in arrays.values())
            for kind in ['animal','day']:
                result=nested(pair,h,kind,arrays,graph,runhash);results.append(result)
                print(json.dumps({'pair':pair,'history':h,'grouping':kind,'animals':len(ids),'status':result['status'],'completed_folds':len(result['folds'])}),flush=True)
    assert len({tuple(r['pair']) for r in results})==10
    (PRIVATE/'all_results_PRIVATE.json').write_text(json.dumps(results,indent=2)+'\n')
    pd.DataFrame(graphrows).to_csv(HERE/'ANATOMICAL_COMPARATOR.csv',index=False)
    receipt['completed_UTC']=datetime.datetime.now(datetime.timezone.utc).isoformat();receipt['status']='COMPLETE';receipt['evaluated_pairs']=10;receipt['directed_relationships']=20
    (PRIVATE/'run_receipt_PRIVATE.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (HERE/'RUN_AGGREGATE.json').write_text(json.dumps({'run_hash':runhash,'status':'COMPLETE','eligible_pairs':10,'directions':20,'histories':[10,5,15],'groupings':['animal','day'],'training_recordings_opened':len(receipt['recording_sources_PRIVATE']),'validation_test_recordings_opened':False,'total_candidate_model_fits':sum(r['fit_count'] for r in results),'rank1_nonconverged_candidate_fits':sum(r['nonconverged_rank1_fits'] for r in results),'maximum_rank1_objective_increase':max(r['max_rank1_objective_increase'] for r in results),'full_sheaf_fitted':False},indent=2)+'\n')
if __name__=='__main__':main()
