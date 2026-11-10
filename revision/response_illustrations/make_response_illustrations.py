"""Draw response-only diagnostics from completed or archived results.

No model is fitted and no uncertainty interval is recomputed here. All
inferential intervals are copied from the frozen result tables. Descriptive
plots of local seasons or annual summaries do not imply independent samples.
Run from any directory with Python, pandas, numpy and matplotlib.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / 'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 14,
    'axes.titlesize': 15, 'axes.labelsize': 14,
    'xtick.labelsize': 12, 'ytick.labelsize': 13,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': False, 'savefig.facecolor': 'white',
    'pdf.fonttype': 42, 'svg.fonttype': 'none',
})
BLUE, RED, GREY, GREEN = '#2877a5', '#c65f51', '#7b7b7b', '#397567'
COUNTRIES = ['Argentina', 'Brazil']
MODES = ['prepared', 'nested_linear']
INITS = ['reset_fc', 'continuous_kc0.2', 'continuous_kc0.4', 'continuous_kc0.6']
INIT_LABELS = ['Seasonal reset', 'Continuous 0.2', 'Continuous 0.4', 'Continuous 0.6']
records, sources = [], set()

def read(relative):
    sources.add(relative)
    return pd.read_csv(ROOT / relative)

def title(ax, s):
    ax.set_title(s, loc='left', fontweight='bold', pad=12)

def save(fig, key, caption, relevant_comments):
    fig.savefig(OUT / f'{key}.png', dpi=300, bbox_inches='tight', pad_inches=.12)
    fig.savefig(OUT / f'{key}.pdf', dpi=300, bbox_inches='tight', pad_inches=.12)
    plt.close(fig)
    records.append(dict(reference=key.replace('_', ' ', 1),
        source=f'response_illustrations/figures/{key}.png',
        caption=caption, relevant_comments=relevant_comments))

def interval(ax, row, y, color='black', marker='o', ms=6):
    ax.errorbar(row.skill_percent, y,
        xerr=[[row.skill_percent-row.lower], [row.upper-row.skill_percent]],
        fmt=marker, color=color, markersize=ms, capsize=3, linewidth=1.5)

# R1: all lifecycle memberships, together with the national conditional test.
C = read('results_lifecycle/lifecycle_counts.csv')
L = read('results_lifecycle/lifecycle_summary.csv')
assert C.harvest_years.sum() == 50
fig, axs = plt.subplots(2, 1, figsize=(10.8, 8.8), gridspec_kw={'height_ratios':[1.45, 1]})
labels = ['El Niño: developing', 'El Niño: persistent', 'El Niño: reversal',
          'La Niña: developing', 'La Niña: persistent', 'La Niña: reversal',
          'Neutral after event', 'Neutral continuation']
for i, row in C.iterrows():
    axs[0].barh(i, row.harvest_years, color=([RED]*3+[BLUE]*3+[GREY]*2)[i], height=.64)
    axs[0].text(row.harvest_years+.2, i, str(row.harvest_years), va='center', fontsize=13)
axs[0].set(yticks=range(8), yticklabels=labels, ylim=(7.6,-.6), xlim=(0,10.8))
axs[0].set_xlabel('Harvest years in 1975–2024')
axs[0].set_xticks([0,2,4,6,8,10])
title(axs[0], 'a) Lifecycle sample sizes')
for ci, country in enumerate(COUNTRIES):
    for mi, model in enumerate(['ENSO_lag','ENSO_son']):
        row=L[(L.country==country)&(L.model==model)&(L.baseline=='ENSO')&(L.scale=='national')].iloc[0]
        interval(axs[1], row, ci*3+mi, BLUE if ci==0 else RED)
axs[1].axvline(0, color='.55', lw=1, ls='--')
axs[1].set(yticks=[0,1,3,4], yticklabels=['Argentina: prior DJF','Argentina: prior SON',
    'Brazil: prior DJF','Brazil: prior SON'], ylim=(4.7,-.7), xlim=(-10,18))
axs[1].set_xlabel('National error reduction beyond concurrent DJF (%)')
axs[1].set_xticks([-10,0,10])
title(axs[1], 'b) National antecedent ENSO increment')
fig.subplots_adjust(left=.29,right=.97,top=.94,bottom=.1,hspace=.65)
save(fig, 'Figure_R1',
    'Figure R1. ENSO lifecycle support and national-scale antecedent sensitivity. '
    '(a) Exact harvest-year counts from Table S14. Developing means a warm or cold event following a neutral year; persistent events retain the phase, reversals switch warm/cold phase, and neutral years distinguish decay from continuation. '
    'The counts sum to 50 harvests in 1975–2024; only one La Niña harvest follows neutral conditions. '
    '(b) Additional national squared-error reduction from previous-year DJF or preceding SON beyond concurrent DJF, using the matched lifecycle validation. Points and bars are the archived estimates and paired 95% three-year block-bootstrap intervals; every interval includes zero. '
    'The field-scale comparisons are reproduced in Figure S7. Counts are harvest years, not independent ENSO episodes. This figure was prepared for this response from the completed S11 outputs.',
    ['Reviewer 1 Major comment 1'])

# R2: original reset versus continuous storage on all archived unit-seasons.
W=read('results_water/water_initialization_local.csv.gz')
D=read('results_water/water_memory_diagnostics.csv')
P=W.pivot(index=['country','spatial_unit_id','harvest_year'], columns='initialization',
          values='reproductive_requirement_satisfaction')
assert len(P)==85840
fig, axs=plt.subplots(2,2,figsize=(11,8.2))
exceedance=[]
hexagons=[]
for ci,country in enumerate(COUNTRIES):
    z=P.xs(country,level='country')
    old=1-z.reset_fc.to_numpy();new=1-z['continuous_kc0.4'].to_numpy()
    ax=axs[0,ci]
    h=ax.hexbin(old,new,gridsize=50,extent=(0,1,0,1),mincnt=1,
                norm=LogNorm(vmin=1,vmax=10000),cmap='Blues',linewidths=0,rasterized=True)
    hexagons.append(h)
    ax.plot([0,1],[0,1],color=RED,lw=1.1,ls='--')
    ax.set(xlim=(0,1),ylim=(0,1),xlabel='Original reset raw CWS',
           ylabel='Continuous raw CWS',xticks=[0,.5,1],yticks=[0,.5,1])
    title(ax,f'{"ab"[ci]}) {country}')
    dr=D[(D.country==country)&(D.initialization=='continuous_kc0.4')].iloc[0]
    assert len(z)==dr.unit_seasons
    assert abs(np.corrcoef(old,new)[0,1]-dr.CWS_correlation)<1e-12
    ax.text(.04,.95,f'n = {len(z):,} unit–seasons\nr = {dr.CWS_correlation:.6f}\nMaximum |Δ| = {dr.maximum_absolute_ratio_difference:.3f}',
            transform=ax.transAxes,va='top',fontsize=12,bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
    thresholds=np.geomspace(1e-6,.2,200)
    for k,color in zip([.2,.4,.6],[GREY,BLUE,RED]):
        delta=(z.reset_fc-z[f'continuous_kc{k}']).to_numpy()
        assert delta.min()>-1e-10
        percent=np.array([(delta>t).mean()*100 for t in thresholds])
        axs[1,ci].plot(thresholds,percent,label=f'Off-season Kc = {k}',color=color,lw=2)
        exceedance.extend(dict(country=country,off_season_kc=k,threshold=t,percent_above=v)
                          for t,v in zip(thresholds,percent))
    axs[1,ci].set(xscale='log',xlim=(1e-6,.2),ylim=(0,30),
        xlabel='Increase in raw CWS (threshold)',ylabel='Unit–seasons above threshold (%)',
        xticks=[1e-6,1e-4,1e-2,1e-1])
    axs[1,ci].axvline(.01,color='.6',ls=':',lw=1)
    title(axs[1,ci],f'{"cd"[ci]}) Distribution of local changes')
for collection in hexagons:
    collection.set_clim(1,max(float(c.get_array().max()) for c in hexagons))
fig.subplots_adjust(left=.1,right=.87,top=.94,bottom=.15,wspace=.4,hspace=.55)
cax=fig.add_axes([.9,.59,.018,.31]);fig.colorbar(h,cax=cax,label='Unit–seasons per hexagon')
fig.legend(*axs[1,0].get_legend_handles_labels(),loc='lower center',ncol=3,frameon=False,fontsize=12)
pd.DataFrame(exceedance).to_csv(OUT/'Figure_R2_exceedance_values.csv',index=False)
save(fig,'Figure_R2',
    'Figure R2. Original seasonal reset compared with continuous soil-water storage. '
    '(a,b) Raw CWS = 1 minus reproductive water-requirement satisfaction for all 85,840 administrative unit–seasons (10,520 in Argentina; 75,320 in Brazil), using off-season Kc = 0.4. '
    'Hexagon colour shows counts on a logarithmic scale; the dashed line is equality. '
    '(c,d) Percentage of unit–seasons exceeding each increase in raw CWS for off-season Kc = 0.2, 0.4 and 0.6; zero changes are retained in each denominator. '
    'These are descriptive comparisons, with no independence assumption for units sharing a weather grid. '
    'The original reset replay agrees to 6.66 × 10⁻¹⁶. Changing initial storage from 100 to 50 mm changes the Kc = 0.4 satisfaction ratio by at most 5.54 × 10⁻⁷ in Argentina and zero in Brazil after spin-up. '
    'The high correlations and small average changes coexist with larger individual differences; this does not establish that real soil-water memory is negligible. Source: completed S12 water outputs and Table S16.',
    ['Reviewer 1 Major comment 1'])

# R3: all 32 ENSO increments, preserving their archived paired intervals.
M=read('results_climate/matched_model_summary.csv')
Q=M[(M.model=='Climate_ENSO')&(M.baseline=='Climate')].copy()
assert len(Q)==32 and (Q.lower<=0).all() and (Q.upper>=0).all()
fig,axs=plt.subplots(2,2,figsize=(11,9.2),sharex=True,sharey=True)
for si,scale in enumerate(['field','national']):
    for ci,country in enumerate(COUNTRIES):
        ax=axs[si,ci]
        for mi,mode in enumerate(MODES):
            for ki,init in enumerate(INITS):
                row=Q[(Q.country==country)&(Q.scale==scale)&(Q.preprocessing==mode)&(Q.initialization==init)].iloc[0]
                interval(ax,row,ki+(-.14 if mi==0 else .14), BLUE if mi==0 else RED,
                         'o' if mi==0 else 's',5)
        ax.axvline(0,color='.5',lw=1,ls='--')
        ax.set(yticks=range(4),yticklabels=INIT_LABELS,ylim=(3.6,-.6),xlim=(-10,25))
        ax.tick_params(axis='y',labelleft=True)
        title(ax,f'{"abcd"[si*2+ci]}) {country} | {scale}')
        ax.set_xlabel('ENSO error reduction\nbeyond climate (%)')
fig.legend(handles=[Line2D([],[],color=BLUE,marker='o',label='Prepared anomalies'),
    Line2D([],[],color=RED,marker='s',label='Train-only linear trends')],
    loc='lower center',ncol=2,frameon=False,fontsize=13)
fig.subplots_adjust(left=.2,right=.98,top=.94,bottom=.12,hspace=.44,wspace=.62)
Q.to_csv(OUT/'Figure_R3_exact_estimates.csv',index=False)
save(fig,'Figure_R3',
    'Figure R3. ENSO increments across all matched preprocessing, initialization and scoring choices. '
    '(a,b) State/province-field scoring; (c,d) national scoring. Each point directly compares Climate+ENSO with Climate-only using identical observations and validation splits within its design. '
    'Points and bars reproduce the 32 estimates and paired 95% three-year block-bootstrap intervals from Tables S17–S20; all intervals include zero. '
    'Reset means each season starts at field capacity; continuous labels give the off-season Kc. '
    'Prepared and train-only designs differ in trend specification and fitting sample, so their difference is not an isolated estimate of information leakage. '
    'The maximum initialization-induced shift in any matched model/baseline score, across the complete 128-row table, is 0.256 percentage points. '
    'These conditional comparisons use observed seasonal weather and do not measure operational preseason forecast skill.',
    ['Reviewer 1 Major comment 1','Reviewer 1 Major comment 2'])

# R4: values as well as unchanged ranks in the matched climate model.
I=read('results_climate/climate_driver_sensitivity.csv')
names={'mean_temperature_c':'Temperature','mean_vpd_kpa':'VPD',
       'precipitation_total_mm':'Precipitation','heavy_rain_fraction':'Heavy rainfall',
       'hot_dry_fraction':'Hot–dry conditions','hot_wet_fraction':'Hot–wet conditions',
       'CWS':'CWS','crop_water_stress':'CWS','water_stress':'CWS'}
print('Climate indices:',I.index_name.unique().tolist())
fig,axs=plt.subplots(2,2,figsize=(11,9.3))
for mi,mode in enumerate(MODES):
    for ci,country in enumerate(COUNTRIES):
        z=I[(I.country==country)&(I.preprocessing==mode)]
        reset=z[z.initialization=='reset_fc'].sort_values('importance_rank')
        ax=axs[mi,ci];ys=np.arange(7)
        for j,(_,r) in enumerate(reset.iterrows()):
            q=z[z['index']==r['index']]
            assert q.importance_rank.nunique()==1
            # Range is an initialization sensitivity span, not a confidence interval.
            ax.plot([q.importance_percent.min(),q.importance_percent.max()],[j,j],
                    color=RED,lw=3,alpha=.8)
            ax.plot(r.importance_percent,j,'o',color='black',mfc='white',ms=8,mew=1.4)
            c=q[q.initialization=='continuous_kc0.4'].iloc[0]
            ax.plot(c.importance_percent,j,'x',color=BLUE,ms=8,mew=1.6)
        labels=[names.get(x,x) for x in reset.index_name]
        assert not any('_' in x for x in labels),labels
        ax.set(yticks=ys,yticklabels=labels,ylim=(6.7,-.7))
        ax.set_xlabel('Permutation error increase\n(% of null-model error)')
        ax.axvline(0,color='.8',lw=.7)
        title(ax,f'{"abcd"[mi*2+ci]}) {country} | '+('prepared' if mi==0 else 'train-only'))
fig.legend(handles=[Line2D([],[],color='black',marker='o',mfc='white',ls='',label='Seasonal reset'),
    Line2D([],[],color=BLUE,marker='x',ls='',label='Continuous Kc = 0.4'),
    Line2D([],[],color=RED,lw=3,label='Span across all initializations')],
    loc='lower center',ncol=3,frameon=False,fontsize=11.5)
fig.subplots_adjust(left=.19,right=.98,top=.94,bottom=.14,hspace=.45,wspace=.7)
save(fig,'Figure_R4',
    'Figure R4. Sensitivity of climate-index importance to soil-water initialization. '
    'Rows order the seven indices by their reset-model importance within each country and preprocessing design. '
    'Open circles show the reset values, blue crosses continuous Kc = 0.4, and red horizontal segments the minimum-to-maximum span over reset and continuous Kc = 0.2, 0.4 and 0.6. '
    'Importance is the held-out squared-error increase after the index is permuted, expressed relative to null-model error. The spans are sensitivity ranges, not confidence intervals; short segments can be hidden by coincident markers. '
    'All seven ranks remain unchanged within each design. The prepared and train-only designs can yield different importance values and orders. '
    'These are the matched country-field climate-model diagnostics from the completed S12 outputs, not a refit of the archived hierarchical local driver map in Figure 2a.',
    ['Reviewer 1 Major comment 1','Reviewer 1 Major comment 2'])

# R5: uncertainty of the primary archived score, kept separate from new designs.
F=read('response_illustrations/inputs/ENSO_multivariate_field_country_summary.csv')
fig,ax=plt.subplots(figsize=(10.7,3.3))
for j,country in enumerate(COUNTRIES):
    r=F[F.country==country].iloc[0]
    v=r.spatial_field_cv_r_squared_percent;lo=r.field_r2_ci_lower_percent;hi=r.field_r2_ci_upper_percent
    ax.errorbar(v,j,xerr=[[v-lo],[hi-v]],fmt='o',color=BLUE if j==0 else RED,capsize=4,lw=2,ms=8)
    ax.text(14.4,j,f'{v:.2f}%  [{lo:.2f}, {hi:.2f}]',va='center',fontsize=13)
ax.axvline(0,color='.5',lw=1,ls='--')
ax.set(yticks=[0,1],yticklabels=COUNTRIES,ylim=(1.6,-.6),xlim=(-3,24),xticks=[0,5,10])
ax.set_xlabel('Primary ENSO field error reduction (%)')
title(ax,'Primary ENSO-only evaluation | 1974–2024')
fig.subplots_adjust(left=.15,right=.99,top=.78,bottom=.29)
save(fig,'Figure_R5',
    'Figure R5. Uncertainty in the primary archived ENSO-only state/province-field scores. '
    'Points show held-out squared-error reductions relative to the primary no-ENSO baseline; bars and bracketed values show the original 95% three-year block-bootstrap intervals. '
    'Argentina: 6.83% (−1.41% to 13.08%); Brazil: 0.07% (−0.80% to 0.85%). Both intervals include zero. '
    'The point-estimate contrast is not a confidence interval or hypothesis test for the between-country difference. '
    'These 1974–2024 results are kept separate from the 1985–2024 matched climate comparisons and the 1975–2024 lifecycle sensitivity. Source: archived primary field summary and Table S8.',
    ['Reviewer 1 Additional comment 1'])

# R6: exactly the rounded entries already reported in supplementary Table S6.
T=read('response_illustrations/inputs/Table_S6_as_reported.csv')
fig,axs=plt.subplots(1,3,figsize=(11,6.3),sharey=True)
labels=[]
for _,r in T.iterrows():
    name=r.macroregion.replace(' Brazil','').replace('Southern','South').replace('Central','Central')
    labels.append(name+' | '+('EN' if r.phase=='El Niño' else 'LN'))
for j,(field,label,ttl) in enumerate([
    ('mean_anomaly_all_and_80','Mean anomaly (t ha⁻¹)','a) Mean anomaly'),
    ('gain_all_and_80','Gain area (%)','b) Gain area'),
    ('cancellation_all_and_80','Cancellation (%)','c) Cancellation')]):
    for i,(_,r) in enumerate(T.iterrows()):
        a,b=[float(s.strip().replace('−','-')) for s in r[field].split('/')]
        axs[j].plot([a,b],[i,i],color='.6',lw=1.2)
        axs[j].plot(a,i,'o',mfc='white',mec='black',ms=7,mew=1.4)
        axs[j].plot(b,i,'s',color=BLUE,ms=6)
    axs[j].set(yticks=range(8),yticklabels=labels,ylim=(7.6,-.6),xlabel=label)
    title(axs[j],ttl)
    if j==0:axs[j].axvline(0,color='.7',ls='--',lw=.8)
axs[0].set_xlim(-.25,.28);axs[0].set_xticks([-.2,0,.2])
axs[1].set_xlim(30,80);axs[2].set_xlim(30,65)
fig.legend(handles=[Line2D([],[],color='black',mfc='white',marker='o',ls='',label='All retained years'),
    Line2D([],[],color=BLUE,marker='s',ls='',label='Years with ≥80% reference-area coverage')],
    loc='lower center',ncol=2,frameon=False,fontsize=12)
fig.subplots_adjust(left=.22,right=.98,top=.88,bottom=.19,wspace=.32)
save(fig,'Figure_R6',
    'Figure R6. Descriptive sensitivity to historical reference-area coverage. '
    'Open circles reproduce all-years estimates and blue squares the ≥80%-coverage subset in Table S6; connectors link the two summaries for the same macroregion and phase. '
    'EN, El Niño; LN, La Niña. The retained/all harvest counts are 14/17 and 16/19 in Argentina, 11/17 and 14/19 in southern Brazil, 10/17 and 12/19 in central Brazil, and 10/16 and 12/17 in MATOPIBA (EN and LN, respectively). '
    'Points reproduce the rounded table entries; no uncertainty intervals are available for this descriptive reaggregation. '
    'The restriction changes the temporal sample. It does not reconstruct historical soybean expansion or isolate a causal coverage effect. The reported mean-anomaly signs remain unchanged.',
    ['Reviewer 1 Additional comment 2'])

# R7: show each actual MATOPIBA harvest, not just the intensity-class mean.
E=read('response_illustrations/inputs/Figure_1_event_spatial_extent.csv')
E=E[E.macroregion.str.contains('MATOPIBA')].copy()
E['intensity']=E.intensity.replace({'Very strong':'Strong'})
fig,axs=plt.subplots(1,2,figsize=(11,4.9),sharey=True)
numeric=[]
for ci,phase in enumerate(['El Niño','La Niña']):
    ax=axs[ci];color=RED if ci==0 else BLUE
    z=E[E.phase==phase]
    for j,intensity in enumerate(['Weak','Moderate','Strong']):
        q=z[z.intensity==intensity].sort_values('analysis_year')
        assert len(q)>0
        values=q.gain_area_percent.to_numpy()
        assert len(q)==([5,4,7] if ci==0 else [7,4,6])[j]
        offsets=np.linspace(-.16,.16,len(q))
        ax.scatter(j+offsets,values,s=42,color=color,alpha=.75,edgecolors='white',linewidths=.5)
        ax.plot(j,values.mean(),'D',color='black',ms=7)
        ax.text(j,105,f'n = {len(q)}',ha='center',fontsize=12)
        numeric.append(dict(phase=phase,intensity=intensity,n=len(q),mean_gain_area_percent=values.mean()))
    ax.set(xticks=range(3),xticklabels=['Weak','Moderate','Strong'],ylim=(-3,112),
           xlim=(-.5,2.5),xlabel='Episode peak intensity class')
    ax.set_yticks([0,25,50,75,100])
    title(ax,f'{"ab"[ci]}) MATOPIBA | {phase}')
axs[0].set_ylabel('Reference area with yield gain (%)')
fig.legend(handles=[Line2D([],[],marker='o',color=GREY,ls='',label='Individual harvest year'),
    Line2D([],[],marker='D',color='black',ls='',label='Mean of annual shares')],
    loc='lower center',ncol=2,frameon=False,fontsize=12)
fig.subplots_adjust(left=.1,right=.98,top=.88,bottom=.22,wspace=.25)
pd.DataFrame(numeric).to_csv(OUT/'Figure_R7_intensity_means.csv',index=False)
save(fig,'Figure_R7',
    'Figure R7. Harvest-to-harvest spread behind the MATOPIBA intensity summaries. '
    'Coloured points are the archived annual reference-area gain shares; black diamonds are their arithmetic means, giving each harvest equal weight as in Figure 1g. '
    'Points are offset horizontally to remain visible; the offset has no quantitative meaning. Strong includes very strong; classes use the full qualifying episode peak, which may follow the focal harvest. '
    'Counts are harvest years, not independent ENSO events. El Niño weak/moderate/strong groups contain 5/4/7 harvests; La Niña contains 7/4/6. '
    'The plot is descriptive and shows no confidence intervals or test of an intensity trend. Uneven groups, variation among harvests and changing coverage limit a monotonic or causal interpretation. '
    'Source: the original Figure 1 annual spatial summaries and Tables S4–S5.',
    ['Reviewer 2 Comment 8'])

hashes={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sorted(sources)}
for r in records:r['sha256']=hashlib.sha256((ROOT/r['source']).read_bytes()).hexdigest()
(HERE/'figure_manifest.json').write_text(json.dumps(dict(figures=records,source_sha256=hashes,
    no_model_refit=True,uncertainty='Frozen archived intervals only; other plots explicitly descriptive'),
    ensure_ascii=False,indent=2)+'\n')
print(f'Created {len(records)} complete response-only figures; no model refit.')
