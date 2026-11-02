from pathlib import Path
import h5py,numpy as np,json
D=Path('data/worm_behavior_sheaf_readiness_v1');out=[]
for p in sorted((D/'recordings').glob('*-cudata.h5')):
 with h5py.File(p,'r') as f:
  tc=f['timing/timestamp_confocal'][:];tn=f['timing/timestamp_nir'][:];sx=f['behavior/stage_x_nir'][:];sy=f['behavior/stage_y_nir'][:];xc=f['behavior/stage_x'][:];yc=f['behavior/stage_y'][:]
  offset=tn[0];T=tn-offset;ct=np.r_[0,np.cumsum(T)];cx=np.r_[0,np.cumsum(sx)];cy=np.r_[0,np.cumsum(sy)];matches=[]
  for k,t in enumerate(tc):
   mid=np.searchsorted(tn,t);found=[]
   for start in range(max(0,mid-16),min(len(tn),mid+1)):
    for n in range(1,17):
     stop=start+n
     if stop>len(tn) or not start<=mid<=stop:continue
     if abs((ct[stop]-ct[start])/n-(t-offset))>2e-7:continue
     if abs((cx[stop]-cx[start])/n-xc[k])<2e-5 and abs((cy[stop]-cy[start])/n-yc[k])<2e-5:found.append([start,stop])
   matches.append(found)
  stats={'animal':p.name.removesuffix('-cudata.h5'),'frames':len(tc),'unique_joint_bins':sum(len(x)==1 for x in matches),'ambiguous_bins':sum(len(x)>1 for x in matches),'unmatched_bins':sum(not x for x in matches),'matched_time_and_stage_means':True}
  out.append({'summary':stats,'matching_bins_zero_based_PRIVATE':matches});print(stats)
(D/'sample_ingestion_v2/reconstructed_clock_bins_PRIVATE.json').write_text(json.dumps(out)+'\n')
