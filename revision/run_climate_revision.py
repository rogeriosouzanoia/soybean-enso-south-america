"""Matched ENSO / Climate / ENSO+Climate comparisons and water sensitivity.

Usage: python run_climate_revision.py --project PATH
Runs a conditional prepared-anomaly comparison and a fold-nested linear-trend
sensitivity. Seven raw climate variables are required for the nested sensitivity.
"""
from pathlib import Path
import argparse, hashlib, json, platform, sys
import numpy as np
import pandas as pd
from revision_models import validate, contiguous_episode_folds, plot_summary

METRICS=['mean_temperature_c','mean_vpd_kpa','precipitation_total_mm',
         'heavy_rain_fraction','hot_dry_fraction','hot_wet_fraction']


def collect_importance(output):
    rows=[]
    for country in ['Argentina','Brazil']:
        for mode in ['prepared','nested_linear']:
            for variant in ['reset_fc','continuous_kc0.4','continuous_kc0.2','continuous_kc0.6']:
                name=f'{country}_{mode}_{variant}'
                z=pd.read_csv(output/f'{name}_climate_importance.csv')
                z['country']=country;z['preprocessing']=mode;z['initialization']=variant
                z['index_name']=z['index'].map(dict(enumerate(METRICS+['CWS'])))
                z['importance_rank']=z.importance_percent.rank(ascending=False,method='min').astype(int)
                rows.append(z)
    result=pd.concat(rows,ignore_index=True)
    result.to_csv(output/'climate_driver_sensitivity.csv',index=False)
    return result


def ids(data):
    d=data.copy()
    value=d.spatial_unit_id.astype(str).str.replace(r'\.0$','',regex=True).str.replace(r'[^0-9]','',regex=True)
    d['unit']=[str(x).zfill(7 if c=='Brazil' else 5) for x,c in zip(value,d.country)]
    return d


def transform(data,groups,column,training,relative=False,standardize=False):
    """Train-only linear trends with >=8 records per group; no imputation."""
    d=data.copy();tr=d[d.year.isin(training)&np.isfinite(d[column])].copy()
    tr['xx']=tr.year-2000.;tr['xy']=tr.xx*tr[column];tr['xsq']=tr.xx**2
    s=tr.groupby(groups).agg(n=(column,'size'),sy=(column,'sum'),sx=('xx','sum'),
                            sxy=('xy','sum'),sxx=('xsq','sum'))
    denom=s.sxx-s.sx*s.sx/s.n
    s['b']=(s.sxy-s.sx*s.sy/s.n)/denom;s['a']=s.sy/s.n-s.b*s.sx/s.n
    s.loc[(s.n<8)|(denom<=0),['a','b']]=np.nan
    d=d.merge(s[['a','b']],left_on=groups,right_index=True,how='left',validate='many_to_one')
    trend=d.a+d.b*(d.year-2000);d['value']=d[column]-trend
    if relative:d['value']=np.where(trend>.05,100*d.value/trend,np.nan)
    if standardize:
        sd=d[d.year.isin(training)].groupby(groups).value.std(ddof=1).rename('sd')
        d=d.merge(sd,left_on=groups,right_index=True,how='left',validate='many_to_one')
        d['value']=np.where(d.sd>1e-12,d.value/d.sd,np.nan)
    return d


def main():
    p=argparse.ArgumentParser();p.add_argument('--project',type=Path,required=True)
    p.add_argument('--water',type=Path,default=Path(__file__).parent/'results_water'/'water_initialization_local.csv.gz')
    p.add_argument('--output',type=Path,default=Path(__file__).parent/'results_climate')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    root=a.project/'climate_mechanisms';source=a.project/'data_merged'/'soybean_ENSO_spatial_panel_1974_2024.csv'
    fallback=Path(__file__).parent/'inputs'/'recovered_enso_local_panel.csv.gz'
    if not source.exists():source=fallback
    files=[source,root/'NASA_POWER_grid_manifest.csv',root/'soybean_climate_indices_grid_year.csv.gz',
           root/'soybean_water_balance_primary_AWC100.csv',a.water]
    missing=[str(f) for f in files if not f.exists()]
    if missing:raise FileNotFoundError('Missing inputs: '+'; '.join(missing))
    local=ids(pd.read_csv(source,dtype={'spatial_unit_id':str})).rename(columns={'analysis_year':'year'})
    if local.duplicated(['country','unit','year']).any():raise ValueError('Duplicate agricultural keys')
    mapping=ids(pd.read_csv(files[1],dtype={'spatial_unit_id':str}))
    mapping=mapping[['country','unit','grid_id','calendar_group']].drop_duplicates()
    if mapping.duplicated(['country','unit']).any():raise ValueError('Ambiguous grid/calendar mapping')
    local=local.merge(mapping,on=['country','unit'],validate='many_to_one')
    enso=pd.read_csv(Path(__file__).parent/'inputs'/'enso_by_year.csv')
    phases=dict(zip(enso.harvest_year,enso.oni_harvest_year_phase))
    fullfold=contiguous_episode_folds(range(1985,2025),phases,5,phases)
    climate=pd.read_csv(files[2]).rename(columns={'harvest_year':'year'})
    climate=climate[(climate.climate_window=='reproductive')&climate.climate_metric.isin(METRICS)&climate.year.between(1985,2024)].copy()
    key=['grid_id','country','calendar_group','year','climate_metric']
    if climate.duplicated(key).any():raise ValueError('Duplicate climate keys')
    if 'climate_value' not in climate:raise ValueError('Raw climate_value is required for train-only sensitivity')
    water=ids(pd.read_csv(a.water,dtype={'spatial_unit_id':str})).rename(columns={'harvest_year':'year'})
    variants=['reset_fc','continuous_kc0.4','continuous_kc0.2','continuous_kc0.6']
    water=water[water.initialization.isin(variants)]
    widewater=water.pivot(index=['country','unit','year'],columns='initialization',values='reproductive_requirement_satisfaction')
    if not set(variants)<=set(widewater.columns):raise ValueError('Incomplete water-initialization variants')
    local=local.merge(widewater,on=['country','unit','year'],validate='one_to_one')
    local=local[local.year.between(1985,2024)].copy()
    # Use exactly the original field's eligible parent regions and fixed weights.
    field=pd.read_csv(Path(__file__).parent/'inputs'/'published_field_oof.csv')
    regioncodes=local.unit.str[:2];local['region']=local.country+'__'+regioncodes
    ref=field.groupby('region_key').reference_area_ha.first()
    local=local[local.region.isin(ref.index)].copy();local['reference']=local.region.map(ref)
    valid_climate=climate.pivot(index=['country','grid_id','calendar_group','year'],columns='climate_metric',values='climate_value').reset_index()
    local=local.merge(valid_climate,on=['country','grid_id','calendar_group','year'],validate='many_to_one')
    finite=['yield_anomaly_percent','yield_harvested_t_ha','harvested_area_ha']+METRICS+variants
    keep=np.isfinite(local[finite].to_numpy()).all(axis=1)&(local.harvested_area_ha>0)
    coverage=local.groupby(['country','region','year']).harvested_area_ha.sum().rename('all_area')
    retained=local.loc[keep].groupby(['country','region','year']).harvested_area_ha.sum().rename('common_area')
    cov=pd.concat([coverage,retained],axis=1);cov['retained_fraction']=cov.common_area/cov.all_area
    cov.to_csv(a.output/'common_support_coverage.csv')
    local=local.loc[keep].copy()
    outputs=[];index_outputs=[];allpred=[]
    for country,c in local.groupby('country',sort=True):
        for mode in ['prepared','nested_linear']:
            for variant in variants:
                print('Matched climate:',country,mode,variant,flush=True)
                def provider(training):
                    z=c.copy()
                    if mode=='nested_linear':
                        y=transform(z,['unit'],'yield_harvested_t_ha',training,relative=True)
                        z['y']=y.set_index(['unit','year']).value.reindex(pd.MultiIndex.from_frame(z[['unit','year']])).to_numpy()
                        for i,metric in enumerate(METRICS):
                            q=climate[(climate.country==country)&(climate.climate_metric==metric)]
                            q=transform(q,['grid_id','calendar_group'],'climate_value',training,standardize=True)
                            q=q[['grid_id','calendar_group','year','value']].rename(columns={'value':f'x{i}'})
                            z=z.merge(q,on=['grid_id','calendar_group','year'],validate='many_to_one')
                    else:
                        z['y']=z.yield_anomaly_percent
                        for i,metric in enumerate(METRICS):
                            q=climate[(climate.country==country)&(climate.climate_metric==metric)]
                            q=q[['grid_id','calendar_group','year','climate_anomaly_standardized']].rename(columns={'climate_anomaly_standardized':f'x{i}'})
                            z=z.merge(q,on=['grid_id','calendar_group','year'],validate='many_to_one')
                    # CWS is negative accumulated ETa/ETc; affine country scaling
                    # is absorbed by train-only field standardization in each fit.
                    z['x6']=-z[variant]
                    columns=['y']+[f'x{i}' for i in range(7)]
                    finite=np.isfinite(z[columns].to_numpy()).all(axis=1)
                    z=z.loc[finite].copy()
                    for col in columns:z[col]=z[col]*z.harvested_area_ha
                    aggregate=z.groupby(['region','year'],sort=True).agg(
                        **{col:(col,'sum') for col in columns},area=('harvested_area_ha','sum'),reference=('reference','first'))
                    for col in columns:aggregate[col]/=aggregate.area
                    aggregate=aggregate.reset_index();aggregate['weight']=aggregate.reference
                    aggregate=aggregate.merge(enso[['harvest_year','oni_djf_c']],left_on='year',right_on='harvest_year',validate='many_to_one')
                    aggregate['oni']=aggregate.oni_djf_c
                    return aggregate
                label=f'{country}_{mode}_{variant}'
                s,pred=validate(provider,range(1985,2025),phases,['ENSO','Climate','Climate_ENSO'],a.output,label,fullfold,inner_contiguous=True)
                s['country']=country;s['preprocessing']=mode;s['initialization']=variant;outputs.append(s)
                pred['country']=country;pred['preprocessing']=mode;pred['initialization']=variant;allpred.append(pred)
                if variant=='reset_fc':
                    # ENSO-to-climate predictions use the identical support and folds.
                    for i,metric in enumerate(METRICS+['CWS']):
                        def index_provider(training):
                            d=provider(training).copy();d['y']=d[f'x{i}'];return d
                        si,_=validate(index_provider,range(1985,2025),phases,['ENSO'],a.output,label+f'_index{i}',fullfold,inner_contiguous=True)
                        si['country']=country;si['preprocessing']=mode;si['index']=metric;index_outputs.append(si)
    summaries=pd.concat(outputs,ignore_index=True);summaries.to_csv(a.output/'matched_model_summary.csv',index=False)
    collect_importance(a.output)
    pd.concat(index_outputs,ignore_index=True).to_csv(a.output/'matched_ENSO_to_climate.csv',index=False)
    predictions=pd.concat(allpred,ignore_index=True);predictions.to_csv(a.output/'all_matched_predictions.csv.gz',index=False)
    # Every model is paired within a run; across water variants check identical support.
    for (country,mode),g in predictions.groupby(['country','preprocessing']):
        sets=[set(zip(z.region,z.year)) for _,z in g.groupby('initialization')]
        if not all(s==sets[0] for s in sets):raise AssertionError('Water variants have unequal evaluation support')
    plot_summary(summaries[(summaries.preprocessing=='nested_linear')&(summaries.initialization=='continuous_kc0.4')],a.output/'Matched_climate_skill.png')
    (a.output/'provenance.json').write_text(json.dumps(dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,
        platform=platform.platform(),input_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
        design='1985–2024 common local support, state/province aggregation, regional ridge slopes, five contiguous outer blocks with intact episodes, one-year buffer, four contiguous inner episode blocks, six ridge penalties; paired 1999 three-year moving-block resamples',
        interpretation='Retrospective conditional comparison using observed weather and contemporaneous ONI. Nested linear sensitivity is a new specification, not a rerun of original trend-family selection or an operational forecast.'),indent=2))
    print('Completed. Return results_climate and results_water for the response letter.',flush=True)


if __name__=='__main__':main()
