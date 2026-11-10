"""Daily water-balance initialization sensitivity using existing NASA files.

Matches the archived FAO56/bucket equations. Continuous simulations include a
prespecified off-season coefficient, NOT observations of off-season land cover.
No weather downloads, source-file writes, or silent interpolation are performed.
"""
from pathlib import Path
import argparse, calendar, json, hashlib
import numpy as np
import pandas as pd


def midpoint(year,a,b):
    ya=year-1 if a>=8 else year;yb=year-1 if b>=8 else year
    start=pd.Timestamp(ya,int(a),1);end=pd.Timestamp(yb,int(b),calendar.monthrange(yb,int(b))[1])
    return start+pd.Timedelta(days=(end-start).days//2)


def stages(n):
    # Python/NumPy and R both use round-to-even for these boundaries.
    a=max(1,int(np.rint(n*.15)));b=min(n,max(a+1,int(np.rint(n*.40))))
    c=min(n,max(b+1,int(np.rint(n*.80))));kc=np.full(n,.4)
    if b>a:kc[a:b]=np.linspace(.4,1.15,b-a+1)[1:]
    kc[b:c]=1.15
    if n>c:kc[c:]=np.linspace(1.15,.5,n-c+1)[1:]
    middle=np.zeros(n,bool);middle[b:c]=True
    return kc,middle


def et0(x,latitude):
    esfn=lambda t:.6108*np.exp(17.27*t/(t+237.3))
    t=x.T2M.to_numpy();tx=x.T2M_MAX.to_numpy();tn=x.T2M_MIN.to_numpy()
    ea=esfn(x.T2MDEW.to_numpy());es=(esfn(tx)+esfn(tn))/2
    delta=4098*esfn(t)/(t+237.3)**2;gamma=.000665*x.PS.to_numpy()
    j=x.index.dayofyear.to_numpy();phi=np.deg2rad(latitude)
    dr=1+.033*np.cos(2*np.pi*j/365);dec=.409*np.sin(2*np.pi*j/365-1.39)
    ws=np.arccos(np.clip(-np.tan(phi)*np.tan(dec),-1,1))
    ra=24*60/np.pi*.0820*dr*(ws*np.sin(phi)*np.sin(dec)+np.cos(phi)*np.cos(dec)*np.sin(ws))
    rs=np.maximum(0,x.ALLSKY_SFC_SW_DWN.to_numpy());rso=.75*ra
    rnl=4.903e-9*((tx+273.16)**4+(tn+273.16)**4)/2*(.34-.14*np.sqrt(np.maximum(0,ea)))*(1.35*np.minimum(1,rs/np.maximum(rso,.01))-.35)
    rn=.77*rs-rnl;u=x.WS2M.to_numpy()
    return np.maximum(0,(.408*delta*rn+gamma*(900/(t+273))*u*np.maximum(0,es-ea))/(delta+gamma*(1+.34*u)))


def bucket(rain,demand,capacity=100.,initial=100.):
    storage=np.full(len(rain),np.nan);eta=storage.copy();s=float(initial)
    for i,(p,e) in enumerate(zip(rain,demand)):
        if not(np.isfinite(p) and np.isfinite(e)):continue
        available=min(capacity,s+max(0.,p));depletion=capacity-available
        ks=1. if depletion<=.5*capacity else np.clip(available/(.5*capacity),0,1)
        eta[i]=min(available,e*ks);s=max(0.,available-eta[i]);storage[i]=s
    return storage,eta


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--project',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=Path(__file__).parent/'results_water')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    root=args.project/'climate_mechanisms'
    manifest=pd.read_csv(root/'NASA_POWER_grid_manifest.csv',dtype={'spatial_unit_id':str})
    calendars=pd.read_csv(root/'soybean_spatial_units_with_calendars.csv',dtype={'spatial_unit_id':str})
    columns=['sowing_start_month','sowing_end_month','harvest_start_month','harvest_end_month']
    units=calendars.merge(manifest[['climate_location_id','grid_id','power_latitude']].drop_duplicates(),
                          on='climate_location_id',validate='one_to_one')
    groups=units[['grid_id','power_latitude','calendar_group']+columns].drop_duplicates()
    if groups.duplicated(['grid_id','calendar_group']).any():raise ValueError('Ambiguous calendars')
    results=[];audit=[];input_hashes={}
    for i,g in enumerate(groups.itertuples(index=False),1):
        checkpoint=args.output/f'water_{g.grid_id}_{g.calendar_group}.csv'
        base=root/'NASA_POWER_daily_grid'/f'{g.grid_id}.csv.gz'
        extra=root/'NASA_POWER_daily_ET0_inputs'/f'{g.grid_id}_ET0.csv.gz'
        # Checkpoint reuse requires identical bytes and configuration.
        hashes={str(p.relative_to(args.project)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [base,extra]}
        input_hashes.update(hashes)
        script_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        key=hashlib.sha256(json.dumps([hashes,g._asdict(),script_hash],sort_keys=True).encode()).hexdigest()
        keyfile=checkpoint.with_suffix('.sha256')
        if checkpoint.exists() and keyfile.exists() and keyfile.read_text()==key:
            results.append(pd.read_csv(checkpoint));continue
        print(f'Water memory {i}/{len(groups)} {g.grid_id} {g.calendar_group}',flush=True)
        x=pd.read_csv(base,parse_dates=['date']).set_index('date')
        q=pd.read_csv(extra,parse_dates=['date']).set_index('date')
        if not x.index.is_unique or not q.index.is_unique:raise ValueError('Duplicate daily dates')
        x=x.join(q[['ALLSKY_SFC_SW_DWN','PS']],validate='one_to_one')
        required=['T2M','T2M_MAX','T2M_MIN','T2MDEW','WS2M','PRECTOTCORR','ALLSKY_SFC_SW_DWN','PS']
        x=x.reindex(pd.date_range('1984-01-01','2024-07-31',freq='D'))
        if not np.isfinite(x[required].to_numpy()).all():
            raise ValueError(f'{g.grid_id}: incomplete daily data from 1984; repair cache before simulating continuous storage')
        x['et0']=et0(x,g.power_latitude)
        if not np.isfinite(x.et0).all():raise ValueError('Nonfinite ET0')
        seasons={};active=np.full(len(x),np.nan)
        for year in range(1985,2025):
            start=midpoint(year,g.sowing_start_month,g.sowing_end_month)
            end=midpoint(year,g.harvest_start_month,g.harvest_end_month)
            ids=np.flatnonzero((x.index>=start)&(x.index<=end))
            if len(ids)!=(end-start).days+1:raise ValueError('Incomplete representative crop cycle')
            kc,middle=stages(len(ids));active[ids]=kc;seasons[year]=(ids,kc,middle)
        rain=x.PRECTOTCORR.to_numpy();e=x.et0.to_numpy();continuous={}
        for k in [.2,.4,.6]:
            demand=np.where(np.isfinite(active),active,k)*e
            st,ea=bucket(rain,demand)
            continuous[f'continuous_kc{k:.1f}']=(st,ea,demand)
        demand=np.where(np.isfinite(active),active,.4)*e
        st,ea=bucket(rain,demand,initial=50.)
        continuous['continuous_kc0.4_initial50']=(st,ea,demand)
        rows=[]
        for year,(ids,kc,middle) in seasons.items():
            demand=kc*e[ids];st,ea=bucket(rain[ids],demand)
            candidates={'reset_fc':(st,ea,demand)}
            candidates.update({name:(s[ids],a[ids],d[ids]) for name,(s,a,d) in continuous.items()})
            for name,(s,a,d) in candidates.items():
                ratio=a[middle].sum()/d[middle].sum()
                if not(-1e-12<=ratio<=1+1e-12):raise AssertionError('Water satisfaction outside [0,1]')
                rows.append(dict(grid_id=g.grid_id,calendar_group=g.calendar_group,harvest_year=year,
                    initialization=name,reproductive_requirement_satisfaction=ratio,
                    sowing_storage_mm=100. if name=='reset_fc' else float(continuous[name][0][ids[0]-1]),
                    season_start_date=str(x.index[ids[0]].date()),season_end_date=str(x.index[ids[-1]].date())))
        z=pd.DataFrame(rows);z.to_csv(checkpoint,index=False);keyfile.write_text(key);results.append(z)
    allrows=pd.concat(results,ignore_index=True)
    unit_keys=units[['country','spatial_unit_id','grid_id','calendar_group']].drop_duplicates()
    out=unit_keys.merge(allrows,on=['grid_id','calendar_group'],validate='many_to_many')
    if out.duplicated(['country','spatial_unit_id','harvest_year','initialization']).any():raise AssertionError('Duplicate unit seasons')
    out.to_csv(args.output/'water_initialization_local.csv.gz',index=False)
    # Exact-reset replay against the original local water-balance file is a gate.
    frozen=pd.read_csv(root/'soybean_water_balance_primary_AWC100.csv',dtype={'spatial_unit_id':str})
    out['spatial_unit_id']=out.spatial_unit_id.astype(str).str.replace(r'\.0$','',regex=True)
    frozen['spatial_unit_id']=frozen.spatial_unit_id.astype(str).str.replace(r'\.0$','',regex=True)
    for d in [out,frozen]:
        d['spatial_unit_id']=[s.zfill(7 if c=='Brazil' else 5) for s,c in zip(d.spatial_unit_id,d.country)]
    pair=out[out.initialization=='reset_fc'].merge(frozen,on=['country','spatial_unit_id','harvest_year'],
        suffixes=('_new','_original'),validate='one_to_one')
    difference=np.abs(pair.reproductive_requirement_satisfaction_new-pair.reproductive_requirement_satisfaction_original)
    report=dict(matched_rows=len(pair),new_reset_rows=int((out.initialization=='reset_fc').sum()),
                original_rows=len(frozen),max_absolute_reset_difference=float(difference.max()),
                reset_replay_passed=bool(len(pair)==len(frozen)==(out.initialization=='reset_fc').sum() and difference.max()<1e-8),
                off_season_assumption='Kc 0.4 primary; 0.2 and 0.6 sensitivities; no inference of actual preceding crops',
                spinup='Continuous daily simulation from 1984-01-01; 1985 first harvest; initial storage 100 vs 50 mm test',
                missing_data='Any missing daily meteorology after 1984 aborts; no interpolation',input_sha256=input_hashes)
    (args.output/'water_replay_audit.json').write_text(json.dumps(report,indent=2))
    if not report['reset_replay_passed']:raise AssertionError('Original water-balance replay differs; inspect audit before interpreting memory results')
    pivot=out.pivot(index=['country','spatial_unit_id','harvest_year'],columns='initialization',values='reproductive_requirement_satisfaction')
    diagnostics=[]
    for country,z in pivot.groupby(level='country'):
        for name in z.columns:
            diff=z[name]-z.reset_fc
            diagnostics.append(dict(country=country,initialization=name,unit_seasons=len(z),
                mean_raw_CWS_difference=float(-diff.mean()),maximum_absolute_ratio_difference=float(diff.abs().max()),
                CWS_correlation=float(z[name].corr(z.reset_fc)),
                initial50_max_difference=float((z['continuous_kc0.4']-z['continuous_kc0.4_initial50']).abs().max())))
    pd.DataFrame(diagnostics).to_csv(args.output/'water_memory_diagnostics.csv',index=False)
    print('Water initialization sensitivity completed and original reset replay verified.',flush=True)


if __name__=='__main__':main()
