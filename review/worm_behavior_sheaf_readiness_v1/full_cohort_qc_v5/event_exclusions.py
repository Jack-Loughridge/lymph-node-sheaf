"""Descriptive disjoint failure decomposition of the unchanged band_events rule."""
from pathlib import Path
import sys,json,collections
import numpy as np,h5py
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'export_recovery_v3'));from audit_behaviour_events import band_events

def decompose(t,distance,radius,band,velocity):
 near=distance>=radius-band*10000;candidates=np.flatnonzero(near & ~np.r_[False,near[:-1]]);counts=collections.Counter(candidate_entries=len(candidates))
 for i in candidates:
  prior=(t>=t[i]-6)&(t<t[i]);target=(t>=t[i]+1)&(t<=t[i]+3)
  if t[i]-t[0]<6 or not prior.any():reason='insufficient_preentry_history'
  elif near[prior].any():reason='not_six_seconds_continuously_outside'
  elif t[i]+3>t[-1]:reason='right_censored_three_second_target'
  elif not target.any() or not np.isfinite(velocity[target]).all():reason='missing_or_nonfinite_velocity_target'
  else:reason='retained'
  counts[reason]+=1
 return dict(counts)
def main():
 private=ROOT/'data/worm_behavior_sheaf_readiness_v1/full_cohort_qc_v5';rows=[]
 for p in sorted((private.parent/'recordings').glob('*-cudata.h5')):
  with h5py.File(p,'r') as f:
   t=f['timing/timestamp_confocal'][:];d=f['behavior/dist_from_cuCenter'][:];r=float(f['behavior/circle_radius'][()]);v=f['behavior/velocity'][:];rev=f['behavior/reversal_vec'][:]
   for band in [.1,.2,.3]:
    counts=decompose(t,d,r,band,v)
    assert counts.get('retained',0)==len(band_events(t,d,r,band,v,rev))
    rows.append({'animal':p.name.removesuffix('-cudata.h5'),'band_mm':band,'counts':counts})
 (private/'event_failure_decomposition_PRIVATE.json').write_text(json.dumps(rows,indent=2)+'\n')
 agg=[]
 for band in [.1,.2,.3]:
  c=collections.Counter()
  for row in rows:
   if row['band_mm']==band:c.update(row['counts'])
  assert c['candidate_entries']==sum(v for k,v in c.items() if k!='candidate_entries')
  agg.append({'band_mm':band,**dict(c)})
 (HERE/'EVENT_EXCLUSIONS_AGGREGATE.json').write_text(json.dumps(agg,indent=2)+'\n');print(json.dumps(agg,indent=2))
if __name__=='__main__':main()
