"""Matched regional-field ridge validation for the AFM soybean ENSO revision.

No network calls. Source data and per-fit audit trails are retained. This module
does not overwrite the submitted analyses. See README_EXECUTAR.txt.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

LAMBDAS = np.r_[0.1, 1., 10., 100., 1000., 10000.]
SEED = 20261104


def episode_folds(years, phases, k, full_phases=None):
    """Allocate intact episodes using the full annual chronology, including gaps."""
    full = full_phases or phases
    ep, last_year, last_phase = 0, None, None
    lookup = {}
    for year in sorted(full):
        phase = full[year]
        if year != (last_year or year - 2) + 1 or phase != last_phase:
            ep += 1
        lookup[year] = ep
        last_year, last_phase = year, phase
    records = pd.DataFrame([(y, phases[y], lookup[y]) for y in sorted(years)],
                           columns=['year', 'phase', 'episode'])
    assignments = {}
    for phase, z in records.groupby('phase', sort=True):
        loads = np.zeros(k, int)
        counts = z.groupby('episode').size()
        for episode, n in sorted(counts.items(), key=lambda v: (-v[1], v[0])):
            fold = int(np.argmin(loads)) + 1
            assignments[episode] = fold
            loads[fold-1] += n
    return {int(r.year): assignments[int(r.episode)] for r in records.itertuples()}


def outside_buffer(years, test):
    excluded = {y+d for y in test for d in (-1, 0, 1)}
    return set(years) - excluded


def contiguous_episode_folds(years, phases, k, full_phases=None):
    """Contiguous temporal blocks, moving boundaries to preserve ENSO episodes."""
    full=full_phases or phases
    lookup={};episode=0;previous=None;previous_year=None
    for year in sorted(full):
        if phases.get(year,full[year])!=previous or previous_year is None or year!=previous_year+1:
            episode+=1
        lookup[year]=episode;previous=full[year];previous_year=year
    selected={}
    for year in sorted(years):selected.setdefault(lookup[year],[]).append(year)
    n=len(years);used=0;folds={}
    for episode,yy in selected.items():
        f=min(k-1,int((used+(len(yy)-1)/2)*k/n))+1
        folds.update({year:f for year in yy});used+=len(yy)
    if set(folds.values())!=set(range(1,k+1)):
        raise ValueError('Too few ENSO episodes for requested contiguous folds')
    return folds


def warm_cold(oni):
    x = np.asarray(oni, float)
    return np.column_stack((np.maximum(x, 0), np.maximum(-x, 0)))


def features(train, valid, model):
    """Fit ENSO scales and lag residualization on unique training years only."""
    e = train[['year', 'oni']].drop_duplicates().sort_values('year')
    scale = e.oni.std(ddof=1)
    scale = scale if np.isfinite(scale) and scale > 1e-12 else 1.
    current, current_v = warm_cold(train.oni)/scale, warm_cold(valid.oni)/scale
    if model == 'ENSO':
        return current, current_v
    if model == 'SON':
        return warm_cold(train.son)/scale, warm_cold(valid.son)/scale
    if model in ('ENSO_lag', 'ENSO_son'):
        source = 'lag' if model == 'ENSO_lag' else 'son'
        unique = train[['year', 'oni', source]].drop_duplicates().sort_values('year')
        base = np.column_stack((np.ones(len(unique)), warm_cold(unique.oni)/scale))
        lag = warm_cold(unique[source])/scale
        coef = np.linalg.lstsq(base, lag, rcond=None)[0]
        r = warm_cold(train[source])/scale - np.column_stack((np.ones(len(train)), current))@coef
        rv = warm_cold(valid[source])/scale - np.column_stack((np.ones(len(valid)), current_v))@coef
        # Do not inflate an almost-zero residual: persistence may remove a term.
        return np.column_stack((current, r)), np.column_stack((current_v, rv))
    climate = [f'x{i}' for i in range(7)]
    x, xv = train[climate].to_numpy(), valid[climate].to_numpy()
    mean, sd = x.mean(axis=0), x.std(axis=0, ddof=1)
    sd = np.where(sd > 1e-12, sd, 1.)
    x, xv = (x-mean)/sd, (xv-mean)/sd
    if model == 'Climate':
        return x, xv
    if model == 'Climate_ENSO':
        return np.column_stack((x, current)), np.column_stack((xv, current_v))
    raise ValueError(model)


def ridge_grid(train, valid, model):
    predictions = np.full((len(LAMBDAS), len(valid)), np.nan)
    null = np.full(len(valid), np.nan)
    # All models use the same regional intercepts and regional slope structure.
    xt, xv = features(train, valid, model)
    for region, ids in valid.groupby('region', sort=True).indices.items():
        mask = train.region.to_numpy() == region
        z = train.loc[mask]
        if len(z) < 10:
            continue
        y = z.y.to_numpy()
        X = np.column_stack((np.ones(len(z)), xt[mask]))
        V = np.column_stack((np.ones(len(ids)), xv[ids]))
        penalty = np.eye(X.shape[1]); penalty[0, 0] = 0
        beta = np.linalg.pinv(X.T@X + LAMBDAS[:, None, None]*penalty,
                             rcond=1e-12)@(X.T@y)
        predictions[:, ids] = beta@V.T
        null[ids] = y.mean()
    return predictions, null


def validate(provider, years, phases, models, output, label, outer_folds=None, inner_contiguous=False):
    """provider(train_years) repeats any train-dependent preprocessing."""
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    years = sorted(years)
    folds = outer_folds or episode_folds(years, phases, 5)
    result, tuning, splits, coefficients, importance = [], [], [], [], []
    for fold in sorted({folds[y] for y in years}):
        te = {y for y in years if folds[y] == fold}
        tr = outside_buffer(years, te)
        if not te:
            continue
        inner = (contiguous_episode_folds if inner_contiguous else episode_folds)(tr, phases, 4, phases)
        losses = {m: np.zeros(len(LAMBDAS)) for m in models}
        loss0 = 0.
        for f in sorted(set(inner.values())):
            vi = {y for y in tr if inner[y] == f}
            ti = outside_buffer(tr, vi)
            data = provider(tuple(sorted(ti)))
            a = data[data.year.isin(ti)].reset_index(drop=True)
            b = data[data.year.isin(vi)].reset_index(drop=True)
            fitted = {m: ridge_grid(a, b, m) for m in models}
            keep = np.isfinite(b.y.to_numpy()) & np.isfinite(b.weight.to_numpy())
            for predictions, null in fitted.values():
                keep &= np.all(np.isfinite(predictions), axis=0) & np.isfinite(null)
            if not keep.any():
                raise ValueError(f'No common validation support: {label}/{fold}/{f}')
            b = b.loc[keep].copy()
            w = b.weight / b.groupby('year').weight.transform('sum')
            y = b.y.to_numpy()
            baseline = next(iter(fitted.values()))[1][keep]
            loss0 += np.sum(w*(y-baseline)**2)
            for m, (p, n) in fitted.items():
                losses[m] += np.sum(w.to_numpy()[None, :]*(y[None, :]-p[:, keep])**2, axis=1)
            splits.append(dict(analysis=label, outer=fold, inner=f,
                               train='|'.join(map(str, sorted(ti))),
                               test='|'.join(map(str, sorted(vi))),
                               buffer_ok=not(set(ti)&{y+d for y in vi for d in (-1,0,1)})))
        selected = {}
        for m in models:
            # Tie preference is prespecified and favours stronger regularization.
            eligible = np.flatnonzero(losses[m] <= losses[m].min()+max(1e-12, 1e-10*loss0))
            selected[m] = int(eligible[-1])
            for j, lam in enumerate(LAMBDAS):
                tuning.append(dict(analysis=label, outer=fold, model=m, penalty=lam,
                                   inner_error=losses[m][j], selected=j==selected[m]))
        data = provider(tuple(sorted(tr)))
        a = data[data.year.isin(tr)].reset_index(drop=True)
        b = data[data.year.isin(te)].reset_index(drop=True)
        fitted = {m: ridge_grid(a, b, m) for m in models}
        keep = np.isfinite(b.y.to_numpy()) & np.isfinite(b.weight.to_numpy())
        for p, n in fitted.values():
            keep &= np.all(np.isfinite(p),axis=0)&np.isfinite(n)
        b = b.loc[keep].copy()
        b['fold'] = fold
        b['Null'] = next(iter(fitted.values()))[1][keep]
        for m, (p, n) in fitted.items():
            b[m] = p[selected[m], keep]
        if 'Climate' in models:
            xt,xv=features(a,b,'Climate')
            beta_by_region={}
            lam=float(LAMBDAS[selected['Climate']])
            for region,idx in b.groupby('region',sort=True).indices.items():
                mask=a.region.to_numpy()==region
                X=np.column_stack((np.ones(mask.sum()),xt[mask]));y=a.loc[mask,'y'].to_numpy()
                penalty=np.eye(8);penalty[0,0]=0
                beta=np.linalg.pinv(X.T@X+lam*penalty,rcond=1e-12)@(X.T@y)
                beta_by_region[region]=beta
                for k in range(7):
                    coefficients.append(dict(analysis=label,fold=fold,region=region,index=k,
                                             field_standardized_coefficient=beta[k+1],penalty=lam))
            w=(b.weight/b.groupby('year').weight.transform('sum')).to_numpy()
            observed=b.y.to_numpy();base=b.Climate.to_numpy()
            null_error=float(np.sum(w*(observed-b.Null.to_numpy())**2))
            rng=np.random.default_rng(SEED+int(fold))
            for k in range(7):
                for repetition in range(20):
                    changed=base.copy()
                    for region,idx in b.groupby('region',sort=True).indices.items():
                        perm=rng.permutation(idx)
                        changed[idx]+=(xv[perm,k]-xv[idx,k])*beta_by_region[region][k+1]
                    increase=float(np.sum(w*((observed-changed)**2-(observed-base)**2)))
                    importance.append(dict(analysis=label,fold=fold,index=k,repetition=repetition,
                                           weighted_error_increase=increase,weighted_null_error=null_error))
        result.append(b)
        splits.append(dict(analysis=label,outer=fold,inner=0,
                           train='|'.join(map(str, sorted(tr))),test='|'.join(map(str, sorted(te))),
                           buffer_ok=not(set(tr)&{y+d for y in te for d in (-1,0,1)})))
    predictions = pd.concat(result, ignore_index=True).sort_values(['year','region'])
    if predictions.duplicated(['region','year']).any():
        raise AssertionError('Duplicate held-out predictions')
    if set(predictions.year.unique()) != set(years):
        raise AssertionError('Incomplete outer-year coverage')
    assert all(s['buffer_ok'] for s in splits)
    predictions.to_csv(output/f'{label}_predictions.csv.gz',index=False)
    pd.DataFrame(tuning).to_csv(output/f'{label}_tuning.csv',index=False)
    pd.DataFrame(splits).to_csv(output/f'{label}_splits.csv',index=False)
    if coefficients:
        pd.DataFrame(coefficients).to_csv(output/f'{label}_climate_coefficients.csv',index=False)
        raw=pd.DataFrame(importance);raw.to_csv(output/f'{label}_climate_permutations.csv',index=False)
        avg=raw.groupby(['index','fold'],as_index=False).agg(
            error_increase=('weighted_error_increase','mean'),null_error=('weighted_null_error','first'))
        avg=avg.groupby('index',as_index=False).agg(error_increase=('error_increase','sum'),null_error=('null_error','sum'))
        avg['importance_percent']=100*avg.error_increase/avg.null_error
        avg.to_csv(output/f'{label}_climate_importance.csv',index=False)
    comparisons = [(m,'Null') for m in models]
    if 'ENSO_lag' in models:
        comparisons += [('ENSO_lag','ENSO'), ('ENSO_son','ENSO')]
    if 'Climate_ENSO' in models:
        comparisons += [('Climate_ENSO','Climate')]
    summary, annual = summarize(predictions, comparisons)
    summary.to_csv(output/f'{label}_summary.csv',index=False)
    annual.to_csv(output/f'{label}_annual_errors.csv',index=False)
    return summary, predictions


def summarize(predictions, comparisons, draws=1999):
    annual = []
    for year,z in predictions.groupby('year',sort=True):
        w=z.weight.to_numpy();w=w/w.sum();y=z.y.to_numpy()
        for m,b in comparisons:
            e=y-z[m].to_numpy();e0=y-z[b].to_numpy()
            annual.append(dict(year=int(year),model=m,baseline=b,
                field_error=np.sum(w*e**2),field_baseline=np.sum(w*e0**2),
                national_error=(w@e)**2,national_baseline=(w@e0)**2,
                regions=len(z)))
    annual=pd.DataFrame(annual);summary=[]
    for (m,b),z in annual.groupby(['model','baseline'],sort=True):
        n=len(z);rng=np.random.default_rng(SEED)
        starts=rng.integers(n,size=(draws,int(np.ceil(n/3))))
        ids=((starts[:,:,None]+np.arange(3))%n).reshape(draws,-1)[:,:n]
        for scale in ['field','national']:
            e=z[f'{scale}_error'].to_numpy();e0=z[f'{scale}_baseline'].to_numpy()
            boot=100*(1-e[ids].sum(axis=1)/e0[ids].sum(axis=1))
            lo,hi=np.quantile(boot,[.025,.975])
            summary.append(dict(model=m,baseline=b,scale=scale,
                skill_percent=100*(1-e.sum()/e0.sum()),lower=lo,upper=hi,
                years=n,first_year=z.year.min(),last_year=z.year.max(),
                minimum_regions=z.regions.min(),maximum_regions=z.regions.max()))
    return pd.DataFrame(summary),annual


def plot_summary(summary, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False})
    countries=list(summary.country.unique());fig,axs=plt.subplots(1,len(countries),figsize=(10,4),squeeze=False)
    for ax,country in zip(axs[0],countries):
        z=summary[(summary.country==country)&(summary.scale=='field')].reset_index(drop=True)
        labels=[f'{r.model} vs {r.baseline}' for r in z.itertuples()]
        y=np.arange(len(z))
        for k,r in enumerate(z.itertuples()):
            ax.plot([r.lower,r.upper],[k,k],color='black',lw=1.2)
            ax.plot(r.skill_percent,k,'o',color='black',ms=4)
        ax.axvline(0,color='0.6',ls='--',lw=.8);ax.set_yticks(y,labels);ax.invert_yaxis()
        ax.set_title(country,fontweight='bold');ax.set_xlabel('Held-out squared-error reduction (%)')
    fig.tight_layout();fig.savefig(path,dpi=220,bbox_inches='tight');plt.close(fig)
