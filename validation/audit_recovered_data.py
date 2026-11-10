from pathlib import Path
import pandas as pd
import numpy as np
import json

ROOT=Path(__file__).resolve().parent
I=ROOT/'inputs'; O=ROOT/'results'
checks=[]
def check(label, passed, value):
    checks.append({'check':label,'passed':bool(passed),'value':value})
    if not passed: raise AssertionError(f'{label}: {value}')

d=pd.read_csv(I/'recovered_enso_local_panel.csv.gz',low_memory=False)
check('recovered_panel_rows',len(d)==82259,len(d))
check('recovered_panel_unique_keys',not d.duplicated(['country','spatial_unit_id','analysis_year']).any(),len(d))
d['region_key']=d.country+'__'+d.apply(lambda r:str(int(r.spatial_unit_id)).zfill(7 if r.country=='Brazil' else 5)[:2],axis=1)
d['area_times_anomaly']=d.harvested_area_ha*d.yield_anomaly_percent
g=d.groupby(['country','region_key','analysis_year']).agg(area=('harvested_area_ha','sum'),numerator=('area_times_anomaly','sum')).reset_index()
g['reconstructed_field_anomaly']=g.numerator/g.area
p=pd.read_csv(I/'published_field_oof.csv')
z=p.merge(g,on=['country','region_key','analysis_year'],validate='one_to_one')
err=np.max(np.abs(z.yield_anomaly_percent-z.reconstructed_field_anomaly))
check('recovered_panel_reproduces_field_targets',len(z)==len(p) and err<1e-9,float(err))

# Independent OLS leave-one-year-out reconstruction of Figure 1.
original=pd.read_csv(I/'published_phase_means.csv')
out=[]
for (country,uid),z in d.groupby(['country','spatial_unit_id'],sort=True):
    x=z.analysis_year.to_numpy(dtype=float)-1999.;y=z.yield_harvested_t_ha.to_numpy()
    X=np.column_stack([np.ones(len(x)),x]); inv=np.linalg.inv(X.T@X)
    fit=X@(inv@(X.T@y));h=np.sum((X@inv)*X,axis=1)
    pred=y-(y-fit)/(1-h)
    anomaly=y-pred
    phase=z.oni_harvest_year_phase.to_numpy()
    for label,old_label in [('El Nino','El Niño'),('La Nina','La Niña')]:
        mask=phase==label
        if mask.any():out.append({'country':country,'spatial_unit_id':uid,'phase':old_label,
            'mean_new':np.mean(anomaly[mask]),'median_new':np.median(anomaly[mask]),
            'percent_new':np.mean(100*anomaly[mask]/pred[mask]),'events_new':int(mask.sum())})
new=pd.DataFrame(out)
z=original.merge(new,on=['country','spatial_unit_id','phase'],validate='one_to_one')
check('phase_mean_rows',len(z)==len(original),len(z))
for old,col in [('mean_yield_anomaly_t_ha','mean_new'),('median_yield_anomaly_t_ha','median_new'),('mean_yield_anomaly_percent','percent_new')]:
    err=np.max(np.abs(z[old]-z[col]));check('phase_'+col,err<1e-6,float(err))
check('phase_event_counts',bool((z.events==z.events_new).all()),int((z.events!=z.events_new).sum()))
summary=[]
for (region,phase),g in z.groupby(['macroregion','phase'],sort=True):
    summary.append({'macroregion':region,'phase':phase,'mean_anomaly_t_ha':np.average(g.mean_new,weights=g.reference_area_ha)})
pd.DataFrame(summary).to_csv(O/'recomputed_phase_means_from_raw.csv',index=False)
pd.DataFrame(checks).to_csv(O/'recovered_data_checks.csv',index=False)
print(json.dumps({'passed':sum(x['passed'] for x in checks),'failed':sum(not x['passed'] for x in checks),'checks':checks},indent=2))
