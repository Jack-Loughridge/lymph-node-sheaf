#!/usr/bin/env python3
"""Check every pilot table manifest, joins, feature roles and energy identities."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from prepare import digest,INPUT_GENES,ANNOTATION_GENES


def validate(out):
 config=json.loads((out/'run_config.json').read_text())
 for m in config['tables']:
  path=out/m['file'];df=pd.read_csv(path)
  assert digest(path)==m['sha256'], path
  assert df.shape==(m['rows'],m['columns']), path
  numeric=df.select_dtypes(include=[np.number]).to_numpy()
  assert not np.isinf(numeric).any(), path
 for p in config['code']:assert digest(p['file'])==p['sha256'], p
 assert digest(Path(__file__).with_name('MODEL_SPECIFICATION.txt'))==config['specification_sha256']
 assert digest(config['cache'])==config['cache_sha256']
 read=lambda name:pd.read_csv(out/name)
 joins=read('data_join_audit.csv');annotations=read('annotation_summary.csv')
 assert joins.unmatched_expression_cells.eq(0).all()
 assert (joins.qc_expression_cells==joins.joined_cells).all()
 assert joins.joined_cells.sum()==config['qc_cells']
 assert len(joins)==config['fovs']==4
 assert annotations.groupby('geo').n_cells.sum().to_dict()==joins.set_index('geo').joined_cells.to_dict()
 coverage=read('allocation_summary.csv')
 assert (coverage.exposed_receivers+coverage.zero_exposure_receivers==coverage.n_receivers).all()
 assert (coverage.exposed_receivers<=coverage.receptor_positive_receivers).all()
 assert (coverage.receivers_with_multiple_senders<=coverage.exposed_receivers).all()
 assert (coverage.max_sender_receiver_distance_px<=coverage.radius_px+1e-8).all()
 assert coverage.loc[coverage.exposed_receivers==0,'positive_flux_edges'].eq(0).all()
 weights=read('response_weights.csv');manifest=read('signed_response_gene_manifest.csv')
 assert not weights.gene.isin(set(INPUT_GENES)|set(ANNOTATION_GENES)).any()
 assert not weights.duplicated(['receiver_proxy','programme','gene']).any()
 assert manifest.loc[manifest.source_usable==0,'exclusion_reason'].str.contains('ambiguous_duplicate_sheet').all()
 source=manifest.loc[manifest.exclusion_reason.isna()]
 derived={}
 for (c,f),group in source.groupby(['receiver_proxy','factor']):
  assert group.gene.is_unique
  for row in group.itertuples():
   programme='BMP' if f in ['BMP2','BMP6'] else f
   k=(c,programme,row.gene)
   derived[k]=derived.get(k,0)+np.sign(row.log2fc)/len(group)/(2 if programme=='BMP' else 1)
 derived={k:v for k,v in derived.items() if abs(v)>1e-14}
 actual={(r.receiver_proxy,r.programme,r.gene):r.weight for r in weights.itertuples()}
 assert set(derived)==set(actual)
 assert all(abs(derived[k]-actual[k])<1e-12 for k in derived)
 metrics=read('fixed_allocation_residuals.csv');energy=read('sheaf_energies.csv')
 for _,group in metrics.groupby(['geo','programme','baseline','trial','beta']):
  total=group.loc[group.receiver_proxy=='all'].iloc[0]
  types=group.loc[group.receiver_proxy!='all']
  assert int(total.n_receivers)==int(types.n_receivers.sum())
  assert np.isclose(total.squared_residual_sum,types.squared_residual_sum.sum(),atol=1e-10)
 assert np.allclose(metrics.mse_fixed_allocation,metrics.squared_residual_sum/metrics.n_receivers)
 assert (energy.energy_relaxed<=energy.energy_fixed+1e-7).all()
 assert (energy.stationarity<1e-6).all()
 assert np.allclose(energy.relaxed_energy_per_receiver,energy.energy_relaxed/energy.n_receivers)
 for stage in ['fixed','relaxed']:
  total=energy[[stage+'_'+b+'_energy' for b in ['allocation','aggregation','response']]].sum(axis=1)
  assert np.allclose(total,energy['energy_'+stage if stage=='fixed' else 'energy_relaxed'])
 assert energy.fixed_allocation_energy.max()<1e-15
 assert energy.fixed_aggregation_energy.max()<1e-15
 reference={(r['receiver_proxy'],r['programme']):r for r in json.loads((out/'exposure_normalisation.json').read_text())}
 for r in energy.itertuples():
  ms=metrics.loc[(metrics.geo==r.geo)&(metrics.programme==r.programme)&(metrics.trial==r.trial)&(metrics.beta==r.beta)&(metrics.receiver_proxy!='all')]
  expected=sum(x.squared_residual_sum/(1+(r.beta/reference[x.receiver_proxy,r.programme]['sd'])**2) for x in ms.itertuples())
  assert np.isclose(expected,r.fixed_response_energy,atol=1e-9)
 for _,g in metrics.loc[metrics.receiver_proxy=='all'].groupby(['geo','programme','trial']):
  g=g.set_index('beta')
  # Fixed residual is quadratic in the prespecified gain, an algebra identity.
  f0,f1,f2=float(g.loc[0,'mse_fixed_allocation']),float(g.loc[.25,'mse_fixed_allocation']),float(g.loc[.5,'mse_fixed_allocation'])
  a=(f2-2*f1+f0)/(.125);b=(f1-f0-a*.0625)/.25
  assert np.isclose(g.loc[1,'mse_fixed_allocation'],a+b+f0,atol=1e-10)
 result=read('residual_summary.csv')
 for r in result.itertuples():
  group=metrics.loc[(metrics.programme==r.programme)&(metrics.beta==r.beta)&(metrics.receiver_proxy=='all')]
  obs=group.loc[group.trial==-1]
  assert r.n_receivers==obs.n_receivers.sum()
  assert np.isclose(r.observed_mse,obs.squared_residual_sum.sum()/obs.n_receivers.sum())
  null=group.loc[group.trial>=0].groupby('trial')[['squared_residual_sum','n_receivers']].sum()
  values=null.squared_residual_sum/null.n_receivers
  assert len(values)==r.null_draws==config['nulls']
  assert np.isclose(r.null_mean_mse,values.mean())
  assert np.isclose(r.null_q025,values.quantile(.025)) and np.isclose(r.null_q975,values.quantile(.975))
 assert np.allclose(result.loc[result.beta==0,'observed_mse'],1,atol=1e-10)
 movement=read('shuffle_movement.csv')
 assert len(movement)==config['fovs']*config['nulls']
 assert movement.fraction_sender_positions_moved.between(0,1).all()
 assert movement.fraction_sender_positions_moved.min()>.8
 report={'table_manifest_checks':'passed','algebra_energy_and_sum_checks':'passed','signed_feature_role_checks':'passed',
  'join_and_annotation_count_checks':'passed','baseline_and_null_recomputation':'passed',
  'source_code_and_cache_hashes':'passed','tables':len(config['tables']),
  'qc_cells':config['qc_cells'],'worst_stationarity':float(energy.stationarity.max()),
  'minimum_shuffle_moved_fraction':float(movement.fraction_sender_positions_moved.min())}
 (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('out',type=Path)
 validate(p.parse_args().out)
