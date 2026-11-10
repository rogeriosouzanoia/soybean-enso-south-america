from pathlib import Path
import json,sys,hashlib
import numpy as np
import pandas as pd

import argparse
parser=argparse.ArgumentParser(description='Audit the returned AFM model outputs without rerunning meteorology.')
parser.add_argument('--input-zip',type=Path,default=None)
args=parser.parse_args()
ROOT=Path(__file__).resolve().parent
B=ROOT
O=ROOT/'results_integration'
O.mkdir(exist_ok=True)
summary=pd.read_csv(B/'results_climate/matched_model_summary.csv')
pred=pd.read_csv(B/'results_climate/all_matched_predictions.csv.gz')
checks=[]
def check(name,condition):
    assert condition,name
    checks.append(name)

check('No duplicate matched held-out keys',not pred.duplicated(['country','preprocessing','initialization','region','year']).any())
check('All held-out data and model predictions finite',np.isfinite(pred[['y','weight','Null','ENSO','Climate','Climate_ENSO']]).all().all())
check('All held-out weights positive',pred.weight.gt(0).all())
recomputed=[]
for (country,mode,initialization),g in pred.groupby(['country','preprocessing','initialization']):
    label=f'{country}_{mode}_{initialization}'
    check(label+' complete forty years',set(g.year)==set(range(1985,2025)))
    splits=pd.read_csv(B/f'results_climate/{label}_splits.csv')
    for row in splits.itertuples():
        tr={int(x) for x in row.train.split('|')};te={int(x) for x in row.test.split('|')}
        check(label+f' buffered split {row.outer}/{row.inner}',not tr&{x+d for x in te for d in [-1,0,1]})
    for m,b in [('ENSO','Null'),('Climate','Null'),('Climate_ENSO','Null'),('Climate_ENSO','Climate')]:
        e=[];e0=[];ne=[];ne0=[]
        for year,z in g.groupby('year',sort=True):
            w=z.weight.to_numpy();w=w/w.sum();err=z.y.to_numpy()-z[m].to_numpy();base=z.y.to_numpy()-z[b].to_numpy()
            e.append(w@(err**2));e0.append(w@(base**2));ne.append((w@err)**2);ne0.append((w@base)**2)
        rng=np.random.default_rng(20261104);starts=rng.integers(40,size=(1999,14))
        ids=((starts[:,:,None]+np.arange(3))%40).reshape(1999,-1)[:,:40]
        for scale,x,y in [('field',e,e0),('national',ne,ne0)]:
            x=np.array(x);y=np.array(y);point=100*(1-x.sum()/y.sum())
            ci=np.quantile(100*(1-x[ids].sum(1)/y[ids].sum(1)),[.025,.975])
            r=summary[(summary.country==country)&(summary.preprocessing==mode)&(summary.initialization==initialization)&(summary.model==m)&(summary.baseline==b)&(summary.scale==scale)].iloc[0]
            check(label+f' {m}/{b} {scale} score and interval',np.allclose([point,*ci],[r.skill_percent,r.lower,r.upper],atol=1e-10,rtol=1e-10))
            recomputed.append(dict(country=country,preprocessing=mode,initialization=initialization,model=m,baseline=b,scale=scale,skill_percent=point,lower=ci[0],upper=ci[1]))
    tune=pd.read_csv(B/f'results_climate/{label}_tuning.csv')
    for key,z in tune.groupby(['outer','model']):
        check(label+f' tuning {key}',z.selected.sum()==1 and z.loc[z.selected,'inner_error'].iloc[0]<=z.inner_error.min()+1e-6)

for (country,mode),g in pred.groupby(['country','preprocessing']):
    sets=[set(zip(z.region,z.year)) for _,z in g.groupby('initialization')]
    check(f'{country}/{mode} identical support over initialization',all(s==sets[0] for s in sets))
    for col in ['y','Null','ENSO']:
        wide=g.pivot(index=['region','year'],columns='initialization',values=col)
        check(f'{country}/{mode} invariant {col}',np.allclose(wide.to_numpy(),wide.iloc[:,0].to_numpy()[:,None],rtol=0,atol=1e-12))

enso=pd.read_csv(ROOT/'inputs/enso_by_year.csv')
phase=dict(zip(enso.harvest_year,enso.oni_harvest_year_phase))
for (country,mode,initialization),g in pred.groupby(['country','preprocessing','initialization']):
    folds=g.groupby('year').fold.first()
    check(f'{country}/{mode}/{initialization} intact outer episodes',all(folds[y]==folds[y-1] for y in range(1986,2025) if phase[y]==phase[y-1]))

water=pd.read_csv(B/'results_water/water_initialization_local.csv.gz',dtype={'spatial_unit_id':str})
check('Water unique unit-season-initialization keys',not water.duplicated(['country','spatial_unit_id','harvest_year','initialization']).any())
check('Water satisfaction bounded',water.reproductive_requirement_satisfaction.between(-1e-12,1+1e-12).all())
check('Water storage bounded',water.sowing_storage_mm.between(-1e-12,100+1e-12).all())
check('Water complete 2146 units by forty years by five variants',len(water)==2146*40*5)
tails=[]
for country,g in water.groupby('country'):
    wide=g.pivot(index=['spatial_unit_id','harvest_year'],columns='initialization',values='reproductive_requirement_satisfaction')
    for init in ['continuous_kc0.2','continuous_kc0.4','continuous_kc0.6']:
        diff=wide.reset_fc-wide[init]
        check(country+'/'+init+' continuous satisfaction not higher than reset',diff.min()>-1e-12)
        tails.append(dict(country=country,initialization=init,n=len(diff),mean_difference=diff.mean(),q50=diff.quantile(.5),q95=diff.quantile(.95),q99=diff.quantile(.99),maximum=diff.max(),fraction_above_001=diff.gt(.01).mean(),fraction_above_005=diff.gt(.05).mean(),correlation=wide[init].corr(wide.reset_fc)))
replay=json.loads((B/'results_water/water_replay_audit.json').read_text())
check('Full reset replay passed',replay['reset_replay_passed'] and replay['matched_rows']==85840 and replay['max_absolute_reset_difference']<1e-8)
pre=json.loads((B/'results_water/preflight/preflight_summary.json').read_text())
check('Preflight complete and no unresolved gaps',pre['passed'] and pre['grid_count']==882 and pre['file_count']==1764 and pre['unresolved_records']==0)
files=pd.read_csv(B/'results_water/preflight/file_audit.csv')
check('All 1764 daily sources audited',len(files)==1764 and files.status.eq('audited').all() and files.unresolved.eq(0).all())
hashes={Path(k.replace('\\','/')).name:v for k,v in replay['input_sha256'].items()}
check('All simulated input hashes match preflight final hashes',all(hashes[Path(s.replace('\\','/')).name]==h for s,h in zip(files.source_file,files.sha256_after)))
rep=pd.read_csv(B/'results_water/preflight/radiation_repair_history.csv')
check('Three recorded automatic off-season radiation estimates',len(rep)==3 and rep.variable.eq('ALLSKY_SFC_SW_DWN').all() and (~rep.within_crop).all() and np.isfinite(rep.estimated_value).all())
check('Final water output hash matches climate input',hashlib.sha256((B/'results_water/water_initialization_local.csv.gz').read_bytes()).hexdigest() in json.loads((B/'results_climate/provenance.json').read_text())['input_sha256'].values())
imp=pd.read_csv(B/'results_climate/climate_driver_sensitivity.csv')
for key,g in imp.groupby(['country','preprocessing','index_name']):
    check('Importance rank stable '+str(key),g.importance_rank.nunique()==1)
pivot=summary.pivot(index=['country','preprocessing','scale','model','baseline'],columns='initialization',values='skill_percent')
shift=pivot.sub(pivot.reset_fc,axis=0).abs().max().to_dict()
pd.DataFrame(recomputed).to_csv(O/'verified_matched_scores.csv',index=False)
pd.DataFrame(tails).to_csv(O/'water_difference_distribution.csv',index=False)
report=dict(status='passed',check_count=len(checks),checks=checks,returned_zip_sha256=hashlib.sha256(args.input_zip.read_bytes()).hexdigest() if args.input_zip else None,max_abs_skill_shift_vs_reset_percentage_points=shift,provenance_note='Raw water_replay_audit missing_data describes the numerical simulator. Before simulation the preflight explicitly interpolated three off-season radiation gaps; the prior manual radiation estimate makes four in total. Raw returned files have been retained unchanged.',scope='Audit of returned data, predictions, tuning and simulation outputs. Original daily weather is not present here; no independent full weather rerun or original hierarchical climate/PCA refit was performed.')
(O/'returned_results_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print('Passed',len(checks),'checks. Maximum score shifts:',shift)
print(pd.DataFrame(tails).to_string(index=False))
