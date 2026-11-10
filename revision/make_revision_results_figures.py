"""Update only Figure 2c and draw Figure S8 from the returned matched scores.

The archived Figure 2 is retained in inputs/Figure_2_before_matched.png so the
driver map and alignment panels remain byte-for-byte identical after compositing.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

ROOT=Path(__file__).resolve().parent
S=pd.read_csv(ROOT/'results_climate/matched_model_summary.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,
    'axes.edgecolor':'black','axes.labelcolor':'black','xtick.color':'black',
    'ytick.color':'black','axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':False,'savefig.facecolor':'white'})
PAIRS=[('ENSO','Null'),('Climate','Null'),('Climate_ENSO','Null'),('Climate_ENSO','Climate')]
LABELS=['ENSO vs null','Climate vs null','Climate + ENSO\nvs null','Climate + ENSO\nvs climate']

def panel(ax,country,mode,offset=0,color='black',marker='o'):
    z=S[(S.country==country)&(S.preprocessing==mode)&(S.initialization=='continuous_kc0.4')&(S.scale=='field')]
    for j,(m,b) in enumerate(PAIRS):
        row=z[(z.model==m)&(z.baseline==b)].iloc[0]
        ax.errorbar(row.skill_percent,j+offset,
            xerr=np.array([[row.skill_percent-row.lower],[row.upper-row.skill_percent]]),
            fmt=marker,color=color,markersize=5,linewidth=1.2,capsize=2)
    ax.axvline(0,color='.55',linestyle='--',linewidth=.8)
    ax.set_yticks(range(4),LABELS)
    ax.set_ylim(3.5,-.5);ax.set_xlim(-50,65);ax.set_xticks([-40,-20,0,20,40,60])
    ax.set_xlabel('Held-out squared-error reduction (%)',fontsize=11)

if __name__=='__main__':
    original=Image.open(ROOT/'inputs/Figure_2_before_matched.png').convert('RGB')
    fig=plt.figure(figsize=(6.5,4.02),dpi=144)
    ax=fig.add_axes([.32,.17,.65,.63])
    panel(ax,'Argentina','nested_linear',-.13,'black','o')
    panel(ax,'Brazil','nested_linear',.13,'.5','^')
    fig.text(.01,.95,'c) Matched country-field prediction',fontsize=15,fontweight='bold',va='top')
    ax.legend(handles=[Line2D([],[],marker='o',color='black',label='Argentina'),
                       Line2D([],[],marker='^',color='.5',label='Brazil')],
              loc='lower left',bbox_to_anchor=(-.42,1.01),ncol=2,frameon=False,fontsize=11)
    p=ROOT/'Figures/Figure_2c_matched.png';fig.savefig(p,dpi=144);plt.close(fig)
    new=Image.open(p).convert('RGB').resize((914,579),Image.Resampling.LANCZOS)
    original.paste(new,(935,1090))
    original.save(ROOT/'Figures/Figure_2.png')
    # Assert that all pixels outside the explicitly replaced panel are preserved.
    before=np.array(Image.open(ROOT/'inputs/Figure_2_before_matched.png').convert('RGB'))
    after=np.array(original);mask=np.ones(before.shape[:2],bool);mask[1090:,935:]=False
    assert np.array_equal(before[mask],after[mask])
    fig,axs=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True)
    for row,mode in enumerate(['prepared','nested_linear']):
        for col,country in enumerate(['Argentina','Brazil']):
            ax=axs[row,col];panel(ax,country,mode)
            label='Prepared anomalies' if mode=='prepared' else 'Train-only linear trends'
            ax.set_title(f'{"abcd"[row*2+col]}) {country} | {label}',fontsize=12,fontweight='bold',loc='left')
            ax.tick_params(axis='y',labelleft=True)
            ax.set_xlim(-50,80);ax.set_xticks([-40,-20,0,20,40,60,80])
    fig.subplots_adjust(left=.18,right=.98,bottom=.08,top=.95,hspace=.4,wspace=.65)
    fig.savefig(ROOT/'Figures/Figure_S8_matched.png',dpi=300);plt.close(fig)
