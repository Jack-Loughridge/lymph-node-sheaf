"""Fixed-grid descriptive event QC on author-supplied radial series; no fit.
Coordinates/contacts remain author-derived, not independently reconstructed.
"""
from pathlib import Path
import json
import h5py,numpy as np
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
PRIVATE=ROOT/'data/worm_behavior_sheaf_readiness_v1/export_recovery_v3'
def band_events(t,distance,radius,band_mm,velocity,reversals):
 t=np.asarray(t);near=np.asarray(distance)>=radius-band_mm*10000
 candidates=np.flatnonzero(near & ~np.r_[False,near[:-1]])
 events=[]
 for i in candidates:
  prior=(t>=t[i]-6)&(t<t[i])
  if t[i]-t[0]<6 or not prior.any() or near[prior].any():continue
  if t[i]+3>t[-1]:continue
  target=(t>=t[i]+1)&(t<=t[i]+3)
  if not target.any() or not np.isfinite(velocity[target]).all():continue
  onset=(np.asarray(reversals)==1)&~np.r_[False,np.asarray(reversals[:-1])==1]
  secondary=(t>t[i])&(t<=t[i]+3)
  events.append(dict(frame_zero_based=int(i),timestamp=float(t[i]),velocity_target=float(np.mean(velocity[target])),at_risk_for_new_reversal=bool(reversals[i]==0),new_reversal_within3s=bool(onset[secondary].any())))
 return events

def main():
 rows=[]
 for p in sorted((ROOT/'data/worm_behavior_sheaf_readiness_v1/recordings').glob('*-cudata.h5')):
  with h5py.File(p,'r') as f:
   t=f['timing/timestamp_confocal'][:];d=f['behavior/dist_from_cuCenter'][:];radius=float(f['behavior/circle_radius'][()]);v=f['behavior/velocity'][:];rev=f['behavior/reversal_vec'][:]
   for band in [.1,.2,.3]:
    events=band_events(t,d,radius,band,v,rev)
    rows.append(dict(animal=p.name.removesuffix('-cudata.h5'),band_mm=band,events=events,retained_event_count=len(events),source='Author radial series; general source stage-unit conversion',independent_head_geometry_verified=False))
 (PRIVATE/'behaviour_band_events_PRIVATE.json').write_text(json.dumps(rows,indent=2)+'\n')
 agg=[]
 for band in [.1,.2,.3]:
  selected=[r for r in rows if r['band_mm']==band];events=[e for r in selected for e in r['events']];risk=[e for e in events if e['at_risk_for_new_reversal']]
  agg.append(dict(band_mm=band,sample_animals=len(selected),animals_with_any_event=sum(r['retained_event_count']>0 for r in selected),animals_with_repeated_events=sum(r['retained_event_count']>=2 for r in selected),event_count_range=[min(r['retained_event_count'] for r in selected),max(r['retained_event_count'] for r in selected)],total_events=len(events),secondary_at_risk_events=len(risk),secondary_positive_events=sum(e['new_reversal_within3s'] for e in risk),geometry_independently_verified=False,not_stimulus_TTL=True))
 (HERE/'BEHAVIOUR_EVENT_AGGREGATE.json').write_text(json.dumps(agg,indent=2)+'\n');print(json.dumps(agg,indent=2))
if __name__=='__main__':main()
