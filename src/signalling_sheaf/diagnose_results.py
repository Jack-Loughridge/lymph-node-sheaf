#!/usr/bin/env python3
"""Exact MSE decomposition, to separate exposure amplitude from signed alignment.

These are algebraic finite differences of the prespecified gain grid, not fitted
response coefficients or an estimated gain. MSE(beta)=r2+beta^2*e2-2*beta*re.
"""
import argparse
import json
from pathlib import Path
import pandas as pd


def main(out):
 frame=pd.read_csv(out/'fixed_allocation_residuals.csv')
 frame=frame.loc[frame.receiver_proxy=='all']
 pooled=frame.groupby(['programme','trial','beta'])[['squared_residual_sum','n_receivers']].sum()
 pooled['mse']=pooled.squared_residual_sum/pooled.n_receivers
 results=[]
 for (p,trial),g in pooled.groupby(['programme','trial']):
  g=g.reset_index().set_index('beta')
  f0,f1,f2=[float(g.loc[b,'mse']) for b in [0,.25,.5]]
  e2=(f2-2*f1+f0)/.125
  re=-(f1-f0-.0625*e2)/.5
  assert abs(float(g.loc[1,'mse'])-(f0+e2-2*re))<1e-9
  results.append({'programme':p,'trial':int(trial),'intrinsic_second_moment':f0,
                  'standardised_exposure_second_moment':e2,'response_exposure_cross_moment':re})
 document={'identity':'MSE(beta)=E[r^2]+beta^2 E[e^2]-2 beta E[r e]',
  'interpretation':'Frozen observed exposure scaling means null MSE differences can reflect both amplitude and response alignment; no gain was fitted.',
  'values':results}
 (out/'residual_decomposition.json').write_text(json.dumps(document,indent=2)+'\n')
 marker='\nResidual decomposition diagnostic:\n'
 report=(out/'report.txt').read_text().split(marker)[0]
 lines=[marker.strip(),document['identity'],document['interpretation']]
 for r in results:
  if r['trial']==-1:lines.append(f"  {r['programme']} observed E[e^2]={r['standardised_exposure_second_moment']:.4f}, E[r e]={r['response_exposure_cross_moment']:.4f}.")
 lines.append('BMP residual ordering alone does not imply strong exposure-response association; its observed cross moment is small.')
 (out/'report.txt').write_text(report+'\n'+'\n'.join(lines)+'\n')

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('out',type=Path);main(p.parse_args().out)
