#!/usr/bin/env python3
"""Small fixed-map spatial sender-receiver sheaf pilot; outputs aggregates only."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
from prepare import prepare,signatures,MARKERS,ANNOTATION_GENES,INPUT_GENES,digest
from allocation import allocate,ROUTES
from algebra import operator,solve

PROGRAMMES=('IFNG','TNF','BMP')


def save_csv(out,name,rows,columns=None):
 frame=pd.DataFrame(rows,columns=columns)
 if not len(frame.columns):raise ValueError('Missing schema for '+name)
 frame.to_csv(out/name,index=False)
 again=pd.read_csv(out/name)
 if len(again)!=len(frame) or list(again.columns)!=list(frame.columns):raise ValueError('CSV round-trip failed')
 return {'file':name,'rows':len(frame),'columns':len(frame.columns),'sha256':digest(out/name)}


def score_responses(data,weights):
 """Outcome normalisation uses pooled selected cells within proxy type; no disease fit."""
 counts=data['counts']; genes=list(data['genes']); labels=data['labels']
 h=np.log1p(1e4*counts/data['total'][:,None])
 ix={g:i for i,g in enumerate(genes)}
 scores={p:np.full(len(h),np.nan) for p in PROGRAMMES}
 references=[];gene_refs=[]
 for (c,p),w in weights.items():
  ids=np.flatnonzero(labels==c)
  row={'receiver_proxy':c,'programme':p,'matched_cells':len(ids),'retained_genes':len(w)}
  if len(ids)<5 or len(w)<5:
   references.append(row|{'scorable':False,'reason':'fewer_than_5_receivers_or_genes'});continue
  g=list(w);values=h[ids][:,[ix[x] for x in g]]
  mu=values.mean(axis=0);sd=values.std(axis=0);valid=sd>1e-8
  z=np.zeros_like(values);z[:,valid]=(values[:,valid]-mu[valid])/sd[valid]
  raw=z@np.array([w[x] for x in g]); rawmean=float(raw.mean());rawsd=float(raw.std())
  for k,gene in enumerate(g):
   gene_refs.append({'receiver_proxy':c,'programme':p,'gene':gene,
                    'log_expression_mean':mu[k],'log_expression_sd':sd[k],'constant_gene':bool(not valid[k])})
  if rawsd<=1e-8:
   references.append(row|{'scorable':False,'reason':'constant_programme'});continue
  scores[p][ids]=(raw-rawmean)/rawsd
  references.append(row|{'scorable':True,'reason':'','raw_score_mean':rawmean,'raw_score_sd':rawsd,'variable_genes':int(valid.sum())})
 return h,scores,references,gene_refs


def permute_bundles(supply,labels,rng):
 perm=np.arange(len(supply))
 for c in sorted(set(labels)):
  ids=np.flatnonzero(labels==c);perm[ids]=rng.permutation(ids)
 return supply[perm],float(np.mean(perm!=np.arange(len(supply))))


def exposure_reference(exposure,labels):
 result={}
 for p,v in exposure.items():
  for c in sorted(set(labels)):
   ids=(labels==c)&np.isfinite(v)
   if not ids.any():continue
   sd=float(v[ids].std())
   result[c,p]={'mean':float(v[ids].mean()),'sd':sd if sd>1e-8 else 1.0,
                'observed_sd':sd,'n':int(ids.sum())}
 return result


def evaluate(r,exposure,mu,sd,beta):
 residual=r-beta*(exposure-mu)/sd
 return float(np.mean(residual**2)),residual


def main(args):
 if args.nulls<2:raise ValueError('Use >=2 shuffles')
 if args.allocator=='commot' and not all(importlib.util.find_spec(m) for m in ('commot','anndata')):
  raise RuntimeError('Requested COMMOT unavailable. Install in an isolated environment or explicitly select --allocator kernel.')
 if args.out.exists():raise FileExistsError('Use a fresh output directory')
 data,info=prepare(args.root,args.signatures,args.cache)
 manifest,weights=signatures(args.signatures)
 h,scores,response_refs,gene_refs=score_responses(data,weights)
 args.out.mkdir(parents=True,exist_ok=False)
 out=args.out;tables=[]
 add=lambda name,rows,columns=None:tables.append(save_csv(out,name,rows,columns))
 add('signed_response_gene_manifest.csv',manifest)
 add('response_weights.csv',[{'receiver_proxy':c,'programme':p,'gene':g,'weight':v}
                              for (c,p),w in weights.items() for g,v in w.items()])
 add('response_normalisation.csv',response_refs)
 add('gene_normalisation.csv',gene_refs)
 (out/'annotation_inspection.json').write_text(json.dumps(info['inspections'],indent=2)+'\n')
 add('data_join_audit.csv',info['joins'])
 groups=[];fovinfo=[];annotations=[];coverage=[];exposure={p:np.full(len(h),np.nan) for p in PROGRAMMES}
 for geo in dict.fromkeys(data['geo']):
  ids=np.flatnonzero(data['geo']==geo)
  if len(ids)>2000:raise ValueError('Pilot FOV exceeds 2000 cells')
  areas=data['areas'][ids]
  if not (np.isfinite(areas)&(areas>0)).all():raise ValueError('Invalid segmentation areas')
  radius=args.radius_diameters*float(np.median(2*np.sqrt(areas/np.pi)))
  lab=data['labels'][ids]
  for c,n in sorted(Counter(lab).items()):
   selected=ids[lab==c]
   annotations.append({'geo':geo,'fov':int(data['fov'][ids[0]]),'receiver_proxy':c,'n_cells':n,
       'median_panck':float(np.median(data['panck'][selected])),
       'median_cd45':float(np.median(data['cd45'][selected])),
       'median_cd68':float(np.median(data['cd68'][selected]))})
  fovinfo.append({'geo':geo,'fov':int(data['fov'][ids[0]]),'n_cells':len(ids),'radius_px':radius,
                 'matched_receiver_proxies':int(np.isin(lab,['TA_like','colonocyte_like']).sum())})
  for p in PROGRAMMES:
   rec=np.flatnonzero(np.isfinite(scores[p][ids]))
   if not len(rec):
    coverage.append({'geo':geo,'programme':p,'n_receivers':0,'status':'no_scorable_receiver_proxies'});continue
   allocation=allocate(h[ids],list(data['genes']),data['coords'][ids],rec,p,radius,args.allocator)
   exposure[p][ids[rec]]=allocation['exposure']
   count_senders=np.zeros(len(rec),int)
   for j in range(len(rec)):
    count_senders[j]=len(set(allocation['send'][allocation['rec']==j]))
   coverage.append({'geo':geo,'programme':p,'n_receivers':len(rec),'status':'scored',
     'receptor_positive_receivers':int((allocation['receptor'].max(axis=1)>0).sum()),
     'exposed_receivers':int((allocation['exposure']>0).sum()),'zero_exposure_receivers':int((allocation['exposure']==0).sum()),
     'receivers_with_multiple_senders':int((count_senders>1).sum()),
     'positive_flux_edges':len(allocation['flux']),'allocated_proxy_mass':float(allocation['flux'].sum()),
     'max_sender_receiver_distance_px':float(allocation['distance'].max(initial=0)),'radius_px':radius})
   groups.append(dict(geo=geo,programme=p,ids=ids,rec=rec,radius=radius,allocation=allocation))
  print(f'{geo}: {len(ids)} QC cells; proxy types {dict(Counter(lab))}',flush=True)
 add('fov_selection.csv',fovinfo);add('annotation_summary.csv',annotations);add('allocation_summary.csv',coverage)
 eref=exposure_reference(exposure,data['labels'])
 (out/'exposure_normalisation.json').write_text(json.dumps([dict(receiver_proxy=c,programme=p,**v) for (c,p),v in eref.items()],indent=2)+'\n')
 metrics=[];energies=[];null_moves=[];algebra_checks=[]
 rng=np.random.default_rng(args.seed)
 geneix={g:i for i,g in enumerate(data['genes'])}
 # Shared four-ligand permutation for all programmes, preserving the input bundle.
 perms={}
 input_cols=[geneix[g] for g in ['IFNG','TNF','BMP2','BMP6']]
 for geo in dict.fromkeys(data['geo']):
  ids=np.flatnonzero(data['geo']==geo)
  for trial in range(args.nulls):
   changed,moved=permute_bundles(h[ids][:,input_cols],data['labels'][ids],rng)
   perms[geo,trial]=changed
   null_moves.append({'geo':geo,'trial':trial,'fraction_sender_positions_moved':moved})
 for group in groups:
  geo,p,ids,rec=group['geo'],group['programme'],group['ids'],group['rec']
  receiver_types=data['labels'][ids[rec]]
  r=scores[p][ids[rec]]
  mu=np.array([eref[c,p]['mean'] for c in receiver_types]);sd=np.array([eref[c,p]['sd'] for c in receiver_types])
  for trial in range(-1,args.nulls):
   baseline='observed' if trial==-1 else 'within_proxy_sender_shuffle'
   if trial==-1:a=group['allocation']
   else:
    shuffled=h[ids].copy();shuffled[:,input_cols]=perms[geo,trial]
    a=allocate(shuffled,list(data['genes']),data['coords'][ids],rec,p,group['radius'],args.allocator)
   for beta in (0,.25,.5,1):
    mse,residual=evaluate(r,a['exposure'],mu,sd,beta)
    for c in ['all',*sorted(set(receiver_types))]:
     select=np.ones(len(r),bool) if c=='all' else receiver_types==c
     metrics.append({'geo':geo,'programme':p,'receiver_proxy':c,'baseline':baseline,'trial':trial,
       'beta':beta,'n_receivers':int(select.sum()),'mse_fixed_allocation':float(np.mean(residual[select]**2)),
       'squared_residual_sum':float(np.sum(residual[select]**2)),
       'exposed_receivers':int((a['exposure'][select]>0).sum())})
    if beta==.5 or (beta==0 and trial==-1):
     op=operator(a['supply'],a['send'],a['rec'],a['route'],a['flux'],r,mu,sd,beta)
     energy=solve(op)
     if trial==-1:
      errors=op.D@op.x0-op.h
      max_identity=max(float(np.max(np.abs(errors[op.blocks[b]]),initial=0)) for b in ['allocation','aggregation'])
      if max_identity>1e-8:raise RuntimeError('Allocation/aggregation do not satisfy construction')
      algebra_checks.append({'geo':geo,'programme':p,'beta':beta,'max_fixed_allocation_aggregation_error':max_identity,
                             'stationarity':energy['stationarity']})
     energies.append({'geo':geo,'programme':p,'baseline':baseline,'trial':trial,'beta':beta,
        'n_receivers':len(r),**energy,'relaxed_energy_per_receiver':energy['energy_relaxed']/len(r)})
  print(f'{geo} {p}: compared {args.nulls} sender-position shuffles and beta grid',flush=True)
 add('fixed_allocation_residuals.csv',metrics);add('sheaf_energies.csv',energies)
 add('shuffle_movement.csv',null_moves);add('algebra_checks.csv',algebra_checks)
 allmetrics=pd.DataFrame(metrics);summary_rows=[]
 for p in PROGRAMMES:
  frame=allmetrics.loc[(allmetrics.programme==p)&(allmetrics.receiver_proxy=='all')]
  for beta in (0,.25,.5,1):
   sub=frame.loc[frame.beta==beta]
   obs=sub.loc[sub.trial==-1]
   if obs.empty:continue
   n=int(obs.n_receivers.sum());actual=float(obs.squared_residual_sum.sum()/n)
   null=sub.loc[sub.trial>=0].groupby('trial').agg({'squared_residual_sum':'sum','n_receivers':'sum'})
   values=null.squared_residual_sum/null.n_receivers
   summary_rows.append({'programme':p,'beta':beta,'n_receivers':n,'observed_mse':actual,
    'null_mean_mse':float(values.mean()),'null_q025':float(values.quantile(.025)),
    'null_q975':float(values.quantile(.975)),
    'observed_lower_tail_rank_fraction':float((1+(values<=actual).sum())/(len(values)+1)),
    'null_draws':len(values)})
 add('residual_summary.csv',summary_rows)
 (out/'parameter_provenance.json').write_text(json.dumps({
  'external_experiments':['signed_log2FC_and_FDR','organoid_receiver_sheet','canonical_receptor_gene_composition'],
  'CosMx_inferred':['normalised_transcripts','spatial_geometry','proxy_annotations','relative_allocations','signed_response_alignment','normalisation_statistics'],
  'assumptions':{'marker_detection_count':1,'marker_quorum':2,'min_total_biological_counts':50,
    'minimum_receivers_per_type':5,'minimum_programme_targets':5,'kernel':'gaussian','kernel_bandwidth':'radius/2',
    'no_self_edges':True,'BMP_shared_capacity_gene':'BMPR2','direction_only_weights':True},
  'uncalibrated':['secretion','receptor_protein','affinity','diffusion_radius','active_ligand_processing','tissue_baseline','exposure_response_gain'],
  'beta_grid':[0,.25,.5,1],'nominal_beta':.5,'maps_fitted_to_response':False,
  'response_feature_exclusions':{'input_genes':INPUT_GENES,'annotation_markers':ANNOTATION_GENES},
  'routes':ROUTES,'marker_sets':MARKERS,
  'reference_statistics':'transductive within proxy type across selected FOVs; no disease labels',
  'spatial_null':'four-ligand bundles permuted within FOV and proxy-type strata; receiver inputs/outcomes fixed',
  'energy_comparability':'Use fixed-allocation MSE as primary; eliminated sheaf metric changes with beta and graph degree'
 },indent=2)+'\n')
 versions={m:importlib.metadata.version(m) for m in ['numpy','pandas','scipy']}
 config={'completed_utc':datetime.now(timezone.utc).isoformat(),'allocator':args.allocator,
  'COMMOT_installed':bool(importlib.util.find_spec('commot')),'nulls':args.nulls,'seed':args.seed,
  'radius_diameters':args.radius_diameters,'fovs':len(fovinfo),'qc_cells':int(len(h)),
  'sources':info['sources'],'versions':versions,'cache':str(args.cache),'cache_sha256':digest(args.cache),
  'code':[{'file':str(p),'sha256':digest(p)} for p in sorted(Path(__file__).parent.glob('*.py'))],
  'specification_sha256':digest(Path(__file__).with_name('MODEL_SPECIFICATION.txt')),
  'tables':tables,'interpretation':'Exploratory fixed-map compatibility under provisional annotations and uncalibrated gains; no disease classifier'}
 (out/'run_config.json').write_text(json.dumps(config,indent=2)+'\n')
 write_report(out,config,summary_rows,coverage,response_refs,energies)
 print('Completed pilot: '+str(out),flush=True)


def write_report(out,config,summary_rows,coverage,response_refs,energies):
 lines=['Minimal IFNG/TNF/BMP signalling compatibility pilot',
  f"{config['qc_cells']} QC cells across {config['fovs']} fixed FOVs; {config['nulls']} seeded within-proxy sender-position shuffles.",
  f"Allocation backend: {config['allocator']}; COMMOT installed: {config['COMMOT_installed']}.",
  'No published local cell-type annotations: TA-like/colonocyte-like receiver assignments are marker proxies.',
  'Stem/goblet mappings remain unavailable because the organoid workbook sheets were duplicated.',
  'Protein channel medians are corroboration only, not calibrated annotation gates.',
  '', 'Primary residual: signed response alignment minus beta times centred/scaled relative incoming allocation.',
  'All annotation markers and input ligand/receptor genes are excluded from response features.',
  'BMP2/BMP6 have separate input routes and a shared BMP response coordinate and BMPR2 capacity.',
  'Normalisation uses selected cells within proxy type: transductive, not held-out validation.',
  'beta=0 is the cell-intrinsic type-intercept baseline; other gains are assumptions, not fitted potency.',
  '', 'Matched response maps:']
 for x in response_refs:lines.append(f"  {x['receiver_proxy']} {x['programme']}: {x['matched_cells']} cells, {x['retained_genes']} retained targets, scorable={x['scorable']}")
 lines+=['','Allocation coverage:']
 for x in coverage:
  if x['status']=='scored':lines.append(f"  {x['geo']} {x['programme']}: {x['exposed_receivers']}/{x['n_receivers']} receivers exposed; {x['receivers_with_multiple_senders']} receive multiple senders; {x['positive_flux_edges']} positive route edges.")
 lines+=['','Nominal beta=.5 pooled fixed-allocation MSE (lower is smaller discrepancy):',
         'programme | intrinsic beta=0 | observed beta=.5 | shuffled mean [descriptive 2.5%,97.5%]']
 for p in PROGRAMMES:
  zero=next((x for x in summary_rows if x['programme']==p and x['beta']==0),None)
  nom=next((x for x in summary_rows if x['programme']==p and x['beta']==.5),None)
  if nom:lines.append(f"{p} | {zero['observed_mse']:.4f} | {nom['observed_mse']:.4f} | {nom['null_mean_mse']:.4f} [{nom['null_q025']:.4f},{nom['null_q975']:.4f}]")
 lines+=['','Interpretation:',
  f"At the nominal gain, {sum(x['observed_mse'] < next(z['observed_mse'] for z in summary_rows if z['programme']==x['programme'] and z['beta']==0)-1e-10 for x in summary_rows if x['beta']==.5)} of the three programmes improve the intrinsic baseline.",
  *[f"{x['programme']}: observed discrepancy is {'smaller' if x['observed_mse']<x['null_mean_mse'] else 'larger'} than shuffled mean; descriptive lower-tail rank fraction {x['observed_lower_tail_rank_fraction']:.4f}." for x in summary_rows if x['beta']==.5],
  'This descriptive comparison does not validate exposure-response prediction.',
  'These findings are conditional on the proxy annotations, sparse receptor coverage and arbitrary gain.',
  'No observed/shuffled difference establishes a disease effect. No gain was selected for best performance.',
  'Beta sensitivity is fully exported. The external perturbation directions do not calibrate tissue exposure.',
  'Unexposed receivers remain in the residual denominator. Do not interpret their baseline expression as signalling failure.',
  'Small receptor-positive counts and FOV truncation can dominate this pilot.',
  '', 'Sheaf construction:',
  'Sender/response variables are boundary anchored; route contributions and receiver totals are auxiliary.',
  'Each allocation, sum, and response factor has a kernel stalk; restrictions are coordinate projections.',
  'D contains unit-norm factor rows. Energy ||D*x-h||^2 is exactly lifted incidence energy after eliminating factor coordinates.',
  'All sender contributions sum before the ONE receiver/programme response constraint.',
  'Sparse least squares adjusts auxiliary fluxes only. This is not exposure calibration or inferred physical transport.',
  'The primary discrepancy keeps the response-blind allocation fixed; relaxed energy is a secondary algebraic diagnostic.',
  'Do not compare sheaf energy across beta values as a predictive metric: normalisation changes the factor metric.',
  f"Negative adjusted flux entries across observed nominal solves: {sum(x['negative_adjusted_fluxes'] for x in energies if x['trial']==-1 and x['beta']==.5)}.",
  f"Negative adjusted exposure entries across observed nominal solves: {sum(x['negative_adjusted_exposures'] for x in energies if x['trial']==-1 and x['beta']==.5)}.",
  f"Worst solver stationarity: {max(x['stationarity'] for x in energies):.3g}.",
  '', 'Limitations / next work:',
  'Obtain validated annotations or a reference-transfer model; review proxy thresholds and tissue morphology.',
  'Validate transfer of TA/colonocyte signed response directions to these tissue cells.',
  'Review background/ambient counts, alternate receptor routes, diffusion-range and normalization sensitivity.',
  'Keep allocation independent of response training. Calibrate exposure-response only with external dose-response data or a justified held-out procedure.',
  'Run larger independent specimens before considering disease comparisons; FOVs are not independent donors.',
  '', 'Reproduce: see src/signalling_sheaf/README.txt. Mathematical definition: MODEL_SPECIFICATION.txt.',
  'Input matrices and cell cache stay local under ignored data/; review outputs contain aggregate tables only.',
  'No source data were deleted or modified, no new large data downloaded, no Git commit or push.',
  '', 'Sources:',
  'Capeling dictionary: https://pmc.ncbi.nlm.nih.gov/articles/PMC12891662/',
  'Marker references: https://pmc.ncbi.nlm.nih.gov/articles/PMC7041720/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9046244/',
  'COMMOT API: https://commot.readthedocs.io/en/latest/api/commot.tl.spatial_communication.html']
 (out/'report.txt').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--root',type=Path,default=Path.home()/'Real_Data/spatial_ibd_sheaf/pilot_6k')
 p.add_argument('--signatures',type=Path,default=Path('outputs/organoid_response_audit_corrected/signed_response_genes.csv'))
 p.add_argument('--cache',type=Path,default=Path('data/signalling_sheaf_inputs.npz'))
 p.add_argument('--out',type=Path,required=True)
 p.add_argument('--allocator',choices=['kernel','commot'],default='kernel')
 p.add_argument('--nulls',type=int,default=20)
 p.add_argument('--seed',type=int,default=20261031)
 p.add_argument('--radius-diameters',type=float,default=2.5)
 args=p.parse_args()
 if args.radius_diameters<=0:p.error('radius must be positive')
 main(args)
