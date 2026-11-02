"""Exportable aggregate-only forest plot; no animal IDs or predictions."""
from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent

def main():
 d=pd.read_csv(HERE/'ALL_DIRECTIONS_PRIMARY.csv');p=pd.read_csv(HERE/'POOLED_CONTRASTS.csv')
 fig,axes=plt.subplots(1,2,figsize=(13,8),gridspec_kw={'width_ratios':[1.75,1]},layout='constrained')
 ax=axes[0];y=np.arange(len(d))
 # Show the full interval, including strong negative findings; no clipped favourable selection.
 ax.errorbar(d.conditional_gain_percent,y,xerr=np.vstack([d.conditional_gain_percent-d.conditional_CI_low,d.conditional_CI_high-d.conditional_gain_percent]),fmt='o',color='#285b7e',capsize=3,markersize=4)
 ax.set_yticks(y,[f'{r.sender} → {r.receiver}  (n={r.independent_training_animals})' for r in d.itertuples()]);ax.invert_yaxis();ax.axvline(0,color='grey',lw=1);ax.axvline(5,color='#b88b37',ls='--',lw=1)
 ax.set_xlabel('Conditional held-out MSE reduction (%)\nPositive favours sender history; 95% animal bootstrap intervals')
 ax.set_title('All 20 directions: beyond receiver + behaviour')
 ax.grid(axis='x',alpha=.2)
 ax=axes[1]
 for shift,contrast,label,color in [(-.08,'VAR_vs_AR','Before behaviour','#286f9b'),(.08,'VAR_beh_vs_AR_beh','After behaviour','#ad4b38')]:
  s=p[(p.grouping=='animal')&(p.comparison==contrast)].sort_values('history');yy=np.arange(3)+shift
  ax.errorbar(s.relative_MSE_reduction_percent,yy,xerr=np.vstack([s.relative_MSE_reduction_percent-s.CI95_low_percent,s.CI95_high_percent-s.relative_MSE_reduction_percent]),fmt='o',color=color,label=label,capsize=4)
 ax.set_yticks(range(3),['5 samples','10 samples (primary)','15 samples']);ax.invert_yaxis();ax.axvline(0,color='grey',lw=1);ax.axvline(5,color='#b88b37',ls='--',lw=1);ax.set_xlim(-1.6,5.4)
 ax.set_xlabel('Pooled MSE reduction (%)\nEqual weights for 14 animals');ax.set_title('Fixed history-length sensitivity');ax.legend(loc='lower right');ax.grid(axis='x',alpha=.2)
 fig.suptitle('Retrospective processed-fluorescence forecasts — training cohort only',fontsize=13)
 fig.savefig(HERE/'FORECAST_DIAGNOSTIC.png',dpi=180);fig.savefig(HERE/'FORECAST_DIAGNOSTIC.pdf');plt.close(fig)
if __name__=='__main__':main()
