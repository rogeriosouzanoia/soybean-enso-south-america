"""Run independent ENSO lifecycle sensitivity from the supplied frozen inputs."""
from pathlib import Path
import argparse, hashlib, json, platform, sys
import numpy as np
import pandas as pd
from revision_models import validate, plot_summary


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--inputs',type=Path,default=Path(__file__).parent/'inputs')
    parser.add_argument('--output',type=Path,default=Path(__file__).parent/'results_lifecycle')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    enso=pd.read_csv(args.inputs/'enso_by_year.csv').sort_values('harvest_year')
    enso['lag']=enso.oni_djf_c.shift(1)
    enso['previous_phase']=enso.oni_harvest_year_phase.shift(1)
    enso['stage']=np.select([
        enso.previous_phase.isna(),
        (enso.oni_harvest_year_phase=='Neutral')&(enso.previous_phase=='Neutral'),
        (enso.oni_harvest_year_phase=='Neutral')&(enso.previous_phase!='Neutral'),
        (enso.oni_harvest_year_phase!='Neutral')&(enso.previous_phase==enso.oni_harvest_year_phase),
        (enso.oni_harvest_year_phase!='Neutral')&(enso.previous_phase=='Neutral')],
        ['Unavailable','Neutral continuation','Decay to neutral','Persistent event','Developing event'],
        default='Phase reversal')
    # These labels use only current and previous DJF: no future phase information.
    enso.to_csv(args.output/'lifecycle_year_membership.csv',index=False)
    counts=enso[enso.harvest_year>=1975].groupby(['oni_harvest_year_phase','stage']).agg(
        harvest_years=('harvest_year','size'),years=('harvest_year',lambda x:'|'.join(map(str,x))))
    counts.to_csv(args.output/'lifecycle_counts.csv')
    fields=pd.read_csv(args.inputs/'published_field_oof.csv')
    fields=fields.rename(columns={'analysis_year':'year','region_key':'region',
                                 'yield_anomaly_percent':'y','score_weight':'weight'})
    fields=fields.merge(enso[['harvest_year','lag','previous_phase','stage','oni_son_previous_c']],
                        left_on='year',right_on='harvest_year',validate='many_to_one')
    fields=fields[fields.year>=1975].copy()
    fields['oni']=fields.oni_djf_c;fields['son']=fields.oni_son_previous_c
    full_phases=dict(zip(enso.harvest_year,enso.oni_harvest_year_phase))
    summaries=[];annuals=[];life_yields=[]
    for country,z in fields.groupby('country',sort=True):
        print('Lifecycle models:',country,flush=True)
        z=z.reset_index(drop=True)
        folds=dict(z[['year','outer_fold']].drop_duplicates().to_numpy())
        summary,pred=validate(lambda train:z,range(1975,2025),full_phases,
            ['ENSO','SON','ENSO_lag','ENSO_son'],args.output,country,folds)
        summary['country']=country;summaries.append(summary)
        # Descriptive stage means are country field aggregates, not independent local replicates.
        rows=[]
        for year,a in z.groupby('year',sort=True):
            w=a.weight/a.weight.sum()
            rows.append(dict(country=country,year=year,yield_anomaly_percent=float(w@a.y),
                             phase=a.official_enso_phase.iloc[0],stage=a.stage.iloc[0],
                             oni=a.oni.iloc[0],lag=a.lag.iloc[0]))
        annual=pd.DataFrame(rows);annuals.append(annual)
        for (phase,stage),a in annual.groupby(['phase','stage'],sort=True):
            life_yields.append(dict(country=country,phase=phase,stage=stage,harvest_years=len(a),
                 mean_prepared_anomaly_percent=a.yield_anomaly_percent.mean(),
                 minimum=a.yield_anomaly_percent.min(),maximum=a.yield_anomaly_percent.max(),
                 years='|'.join(map(str,a.year))))
    summary=pd.concat(summaries,ignore_index=True)
    summary.to_csv(args.output/'lifecycle_summary.csv',index=False)
    pd.concat(annuals).to_csv(args.output/'lifecycle_annual_yields.csv',index=False)
    pd.DataFrame(life_yields).to_csv(args.output/'lifecycle_descriptive_means.csv',index=False)
    plot_summary(summary,args.output/'Figure_S7_lifecycle.png')
    correlation=float(enso.loc[enso.harvest_year>=1975,['oni_djf_c','lag']].corr().iloc[0,1])
    manifest=dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,platform=platform.platform(),
        current_previous_DJF_correlation=correlation,
        analysis='1975–2024; matched frozen field targets; asymmetric ENSO terms; six ridge penalties; original outer episode folds; four inner episode folds; one-year buffers; train-only lag residualization; 1999 paired three-year circular block resamples',
        limitation='Conditional on prepared targets and fixed reference footprint. Added lag terms assess information beyond persistence, not a causal memory effect or operational lead-time forecast.',
        input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.inputs.glob('*.csv')})
    (args.output/'provenance.json').write_text(json.dumps(manifest,indent=2))
    print(summary.to_string(index=False),flush=True)


if __name__=='__main__':main()
