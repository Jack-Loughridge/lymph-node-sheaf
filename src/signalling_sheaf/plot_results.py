#!/usr/bin/env python3
"""Export aggregate gain-sensitivity figure, not a disease comparison."""
import argparse
from pathlib import Path
import os
import json
os.environ.setdefault('MPLCONFIGDIR','/tmp/signalling_sheaf_matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def main(out):
 df=pd.read_csv(out/'residual_summary.csv')
 fig,axes=plt.subplots(1,3,figsize=(11,3.8),sharey=True)
 for ax,p in zip(axes,['IFNG','TNF','BMP']):
  d=df.loc[df.programme==p].sort_values('beta')
  ax.fill_between(d.beta,d.null_q025,d.null_q975,color='#de8f05',alpha=.18,label='Shuffle 2.5–97.5%')
  ax.plot(d.beta,d.null_mean_mse,'s--',color='#de8f05',label='Shuffle mean')
  ax.plot(d.beta,d.observed_mse,'o-',color='#0173b2',label='Observed')
  ax.axhline(1,color='#666666',lw=1,ls=':',label='Intrinsic baseline')
  ax.set_title(p if p!='BMP' else 'BMP2/BMP6 shared programme')
  ax.set_xlabel('Assumed gain β');ax.set_xticks([0,.25,.5,1])
  ax.spines[['top','right']].set_visible(False)
 axes[0].set_ylim(0,1.08*max(df.null_q975.max(),df.observed_mse.max()))
 axes[0].set_ylabel('Fixed-allocation response MSE')
 handles,labels=axes[0].get_legend_handles_labels()
 fig.legend(handles,labels,loc='lower center',ncol=4,frameon=False,bbox_to_anchor=(.5,.05),fontsize=8)
 config=json.loads((out/'run_config.json').read_text())
 n=int(df.n_receivers.max())
 fig.suptitle(f"{config['fovs']} FOV pilot: {n} provisional epithelial receivers",fontsize=12)
 fig.text(.5,.02,f"Gains are uncalibrated assumptions; {config['nulls']} within-FOV/proxy sender shuffles. No disease inference.",ha='center',fontsize=8)
 fig.tight_layout(rect=[0,.17,1,.93])
 fig.savefig(out/'gain_sensitivity.png',dpi=180)
 fig.savefig(out/'gain_sensitivity.pdf')
 plt.close(fig)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('out',type=Path);main(p.parse_args().out)
