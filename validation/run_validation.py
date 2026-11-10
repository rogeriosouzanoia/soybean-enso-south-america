"""Independent ENSO-only validation. Run: python run_validation.py

Inputs and protocol are frozen alongside this script. The linear sensitivities
do not refit the article's original eight-candidate trend-selection workflow.
"""
from pathlib import Path
from functools import lru_cache
import json, hashlib, platform, sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'inputs'
OUTPUT = ROOT / 'results'
OUTPUT.mkdir(exist_ok=True)
YEARS = np.arange(1974, 2025)
LAMBDAS = np.r_[0., 10. ** np.linspace(-4, 4, 41)]
CANDIDATES = [(f, float(l)) for f in ['symmetric', 'asymmetric'] for l in LAMBDAS]
AUDIT = []
CHECKS = []
TUNING = []


def check(name, passed, value=''):
    CHECKS.append({'check': name, 'passed': bool(passed), 'value': value})
    if not passed:
        raise AssertionError(f'{name}: {value}')


def excluded(test):
    return {int(y + d) for y in test for d in [-1, 0, 1]}


def event_folds(years, phases, k):
    years = np.array(sorted(years), dtype=int)
    records = []
    episode = 0
    previous = None
    for y in years:
        phase = phases[int(y)]
        if phase != previous:
            episode += 1
        records.append((int(y), phase, episode))
        previous = phase
    d = pd.DataFrame(records, columns=['analysis_year', 'phase', 'episode'])
    folds = {}
    for phase, z in d.groupby('phase', sort=True):
        lengths = z.groupby('episode').size()
        load = np.zeros(k, dtype=int)
        for ep, length in sorted(lengths.items(), key=lambda v: (-v[1], v[0])):
            f = int(np.argmin(load))
            load[f] += length
            folds[int(ep)] = f + 1
    return {int(r.analysis_year): folds[int(r.episode)] for r in d.itertuples()}


def design(oni, family, scale):
    oni = np.asarray(oni)
    if family == 'symmetric':
        return np.column_stack([np.ones(len(oni)), oni / scale])
    return np.column_stack([np.ones(len(oni)), np.maximum(oni, 0) / scale,
                            np.maximum(-oni, 0) / scale])


def prediction_grid(train, valid):
    """Exactly the archived regional ridge design and SVD cutoff."""
    pred = np.full((len(CANDIDATES), len(valid)), np.nan)
    neutral = pred.copy()
    null = np.full(len(valid), np.nan)
    unique_oni = np.unique(train.oni_djf_c.to_numpy())
    scale = np.std(unique_oni, ddof=1) if len(unique_oni) > 1 else 1.
    if not np.isfinite(scale) or scale <= 0:
        scale = 1.
    for region, ids in valid.groupby('region_key', sort=True).indices.items():
        z = train[train.region_key == region]
        if len(z) < 10:
            continue
        y = z.yield_anomaly_percent.to_numpy()
        null[ids] = np.mean(y)
        for j, family in enumerate(['symmetric', 'asymmetric']):
            X = design(z.oni_djf_c, family, scale)
            Xv = design(valid.iloc[ids].oni_djf_c, family, scale)
            P = np.eye(X.shape[1]); P[0, 0] = 0
            A = X.T @ X + LAMBDAS[:, None, None] * P
            inverse = np.linalg.pinv(A, rcond=X.shape[1] * np.sqrt(np.finfo(float).eps))
            beta = inverse @ (X.T @ y)
            sl = slice(j * len(LAMBDAS), (j + 1) * len(LAMBDAS))
            pred[sl, ids] = beta @ Xv.T
            neutral[sl, ids] = beta[:, 0, None]
    return pred, null, neutral


def tune(splits, label, country, outer_id):
    model_sse = np.zeros(len(CANDIDATES))
    null_sse = 0.
    completed = 0
    for train, valid in splits:
        p, n, _ = prediction_grid(train, valid)
        y = valid.yield_anomaly_percent.to_numpy()
        w = valid.score_weight.to_numpy()
        # Common finite support for every model and baseline.
        ok = np.isfinite(n) & np.all(np.isfinite(p), axis=0) & np.isfinite(y) & (w > 0)
        if not np.any(ok):
            continue
        # Give each represented validation year unit total country weight.
        # This matters when a young region lacks ten training observations.
        w = w.copy()
        yr = valid.analysis_year.to_numpy()
        for year in np.unique(yr[ok]):
            ids = ok & (yr == year)
            w[ids] /= w[ids].sum()
        denominator = np.sum(w[ok] * (y[ok] - n[ok]) ** 2)
        if denominator <= 0:
            continue
        model_sse += np.sum(w[None, ok] * (y[None, ok] - p[:, ok]) ** 2, axis=1)
        null_sse += denominator
        completed += 1
    if completed != len(splits) or null_sse <= 0:
        raise ValueError(f'Incomplete tuning: {label}/{country}/{outer_id}: {completed}')
    score = 1 - model_sse / null_sse
    eligible = np.flatnonzero(score >= score.max() - .0025)
    selected = min(eligible, key=lambda i: (CANDIDATES[i][0] != 'symmetric',
                                           -CANDIDATES[i][1], -score[i]))
    for i, (family, lam) in enumerate(CANDIDATES):
        TUNING.append({'analysis': label, 'country': country, 'outer_id': outer_id,
                       'family': family, 'lambda': lam, 'inner_skill_percent': 100 * score[i],
                       'selected': i == selected, 'inner_splits': completed})
    return int(selected)


def output_prediction(train, valid, selected, label, outer_id):
    p, n, neutral = prediction_grid(train, valid)
    r = valid.copy()
    r['model_prediction_percent'] = p[selected]
    r['null_prediction_percent'] = n
    r['neutral_prediction_percent'] = neutral[selected]
    r['selected_family'], r['selected_lambda'] = CANDIDATES[selected]
    r['analysis'] = label
    r['outer_id'] = outer_id
    return r[np.isfinite(r.model_prediction_percent) & np.isfinite(r.null_prediction_percent)]


def register_split(label, country, outer_id, inner_id, train, test, chronological=False):
    tr, te = sorted(map(int, train)), sorted(map(int, test))
    no_overlap = not set(tr).intersection(excluded(te))
    time_ok = not chronological or max(tr) <= min(te) - 2
    AUDIT.append({'analysis': label, 'country': country, 'outer_id': outer_id,
                  'inner_id': inner_id, 'train_years': '|'.join(map(str, tr)),
                  'test_years': '|'.join(map(str, te)), 'one_year_gap_passed': no_overlap,
                  'chronological_required': chronological, 'chronology_passed': time_ok})
    check(f'split_{label}_{country}_{outer_id}_{inner_id}', no_overlap and time_ok)


class CountryPanel:
    def __init__(self, country, raw, archived):
        self.country = country
        self.archived = archived[archived.country == country].copy()
        self.regions = sorted(self.archived.region_key.unique())
        self.phases = dict(self.archived[['analysis_year', 'official_enso_phase']].drop_duplicates().values)
        self.oni = dict(self.archived[['analysis_year', 'oni_djf_c']].drop_duplicates().values)
        self.folds = dict(self.archived[['analysis_year', 'outer_fold']].drop_duplicates().values)
        self.reference = self.archived.groupby('region_key').reference_area_ha.first().to_dict()
        self.names = self.archived.groupby('region_key').parent_region_name.first().to_dict()
        x = raw[raw.country == country].copy()
        width = 7 if country == 'Brazil' else 5
        x['unit'] = x.spatial_unit_id.astype(int).astype(str).str.zfill(width)
        x['region_key'] = country + '__' + x.unit.str[:2]
        x = x[x.region_key.isin(self.regions)]
        check(f'unique_raw_keys_{country}', not x.duplicated(['unit', 'analysis_year']).any())
        self.units = sorted(x.unit.unique())
        self.unit_regions = x.groupby('unit').region_key.first().reindex(self.units).to_numpy()
        self.Y = x.pivot(index='unit', columns='analysis_year', values='yield_harvested_t_ha').reindex(index=self.units, columns=YEARS).to_numpy()
        self.A = x.pivot(index='unit', columns='analysis_year', values='harvested_area_ha').reindex(index=self.units, columns=YEARS).to_numpy()
        self.valid = np.isfinite(self.Y) & np.isfinite(self.A) & (self.A > 0)
        self.x = YEARS - 1999.
        n, intercept, slope, mx, ssx = self.linear_fit(tuple(YEARS))
        fitted = intercept[:, None] + slope[:, None] * self.x
        h = 1 / n[:, None] + (self.x - mx[:, None]) ** 2 / ssx[:, None]
        loo_trend = self.Y - (self.Y - fitted) / (1 - h)
        self.loo = np.where(self.valid & (loo_trend > .05),
                            100 * (self.Y - loo_trend) / loo_trend, np.nan)

    def linear_fit(self, train_years):
        mask = np.isin(YEARS, train_years)[None, :] & self.valid
        n = mask.sum(axis=1)
        x = self.x[None, :]
        y = np.where(mask, self.Y, 0.)
        sx = (mask * x).sum(axis=1)
        mx = np.divide(sx, n, out=np.zeros_like(sx), where=n > 0)
        my = np.divide(y.sum(axis=1), n, out=np.zeros_like(sx), where=n > 0)
        ssx = (mask * x ** 2).sum(axis=1) - n * mx ** 2
        sxy = (y * x).sum(axis=1) - n * mx * my
        slope = np.divide(sxy, ssx, out=np.zeros_like(sx), where=ssx > 0)
        intercept = my - slope * mx
        return n, intercept, slope, mx, ssx

    @lru_cache(maxsize=256)
    def prepare(self, train_years, mode, forward_weights=False):
        n, intercept, slope, _, _ = self.linear_fit(train_years)
        trend = intercept[:, None] + slope[:, None] * self.x
        mask = self.valid & (n[:, None] >= 8) & (trend > .05)
        if mode == 'linear_loo_control':
            anomalies = self.loo
            mask = self.valid & np.isfinite(anomalies)
        else:
            anomalies = np.divide(100 * (self.Y - trend), trend,
                                  out=np.full_like(self.Y, np.nan), where=trend > .05)
        records = []
        for region in self.regions:
            ids = self.unit_regions == region
            ok = mask[ids] & np.isfinite(anomalies[ids])
            a = np.where(ok, self.A[ids], 0.)
            total = a.sum(axis=0)
            numerator = (a * np.where(ok, anomalies[ids], 0.)).sum(axis=0)
            value = np.divide(numerator, total, out=np.full(len(YEARS), np.nan), where=total > 0)
            if forward_weights:
                last = max(train_years)
                recent = np.isin(YEARS, [y for y in train_years if y >= last - 9]) & (total > 0)
                ref = float(np.mean(total[recent])) if np.any(recent) else np.nan
            else:
                ref = self.reference[region]
            for j, year in enumerate(YEARS):
                if np.isfinite(value[j]) and np.isfinite(ref) and ref > 0:
                    records.append((self.country, region, self.names[region], int(year),
                                    value[j], total[j], int(ok[:, j].sum()), self.oni[int(year)],
                                    self.phases[int(year)], ref))
        d = pd.DataFrame(records, columns=['country', 'region_key', 'parent_region_name',
                         'analysis_year', 'yield_anomaly_percent', 'harvested_area_ha',
                         'local_units', 'oni_djf_c', 'official_enso_phase', 'reference_area_ha'])
        d['score_weight'] = d.reference_area_ha / d.groupby('analysis_year').reference_area_ha.transform('sum')
        return d


def blocked(panel, mode):
    label = mode
    pieces = []
    for fold in sorted(set(panel.folds.values())):
        te = {y for y, f in panel.folds.items() if f == fold}
        tr = set(map(int, YEARS)) - excluded(te)
        register_split(label, panel.country, fold, 0, tr, te)
        outer = panel.archived if mode == 'published_replay' else panel.prepare(tuple(sorted(tr)), mode)
        inner_folds = event_folds(tr, panel.phases, 4)
        splits = []
        for inner in range(1, 5):
            vi = {y for y, f in inner_folds.items() if f == inner}
            ti = tr - excluded(vi)
            register_split(label, panel.country, fold, inner, ti, vi)
            data = panel.archived if mode == 'published_replay' else panel.prepare(tuple(sorted(ti)), mode)
            splits.append((data[data.analysis_year.isin(ti)], data[data.analysis_year.isin(vi)]))
        selected = tune(splits, label, panel.country, fold)
        pieces.append(output_prediction(outer[outer.analysis_year.isin(tr)],
                                       outer[outer.analysis_year.isin(te)], selected, label, fold))
    return pd.concat(pieces, ignore_index=True)


def forward(panel):
    label = 'linear_forward'
    pieces = []
    for year in range(1995, 2025):
        te = {year}; tr = set(range(1974, year - 1))
        register_split(label, panel.country, year, 0, tr, te, True)
        outer = panel.prepare(tuple(sorted(tr)), label, True)
        end = max(tr)
        splits = []
        for inner, start in enumerate([end - 5, end - 3, end - 1], 1):
            vi = {start, start + 1}; ti = {y for y in tr if y <= start - 2}
            register_split(label, panel.country, year, inner, ti, vi, True)
            d = panel.prepare(tuple(sorted(ti)), label, True)
            splits.append((d[d.analysis_year.isin(ti)], d[d.analysis_year.isin(vi)]))
        selected = tune(splits, label, panel.country, year)
        pred = output_prediction(outer[outer.analysis_year.isin(tr)],
                                 outer[outer.analysis_year.isin(te)], selected, label, year)
        pred['maximum_training_year'] = max(tr)
        pieces.append(pred)
    return pd.concat(pieces, ignore_index=True)


def summarize(d):
    annual = []
    for (analysis, country, year), z in d.groupby(['analysis', 'country', 'analysis_year'], sort=True):
        w = z.score_weight.to_numpy(); w = w / w.sum()
        y = z.yield_anomaly_percent.to_numpy()
        p = z.model_prediction_percent.to_numpy(); n = z.null_prediction_percent.to_numpy()
        annual.append({'analysis': analysis, 'country': country, 'analysis_year': int(year),
                       'regions': len(z), 'local_units': int(z.local_units.sum()),
                       'field_model_sse': np.sum(w * (y - p) ** 2),
                       'field_null_sse': np.sum(w * (y - n) ** 2),
                       'observed_national_percent': w @ y,
                       'predicted_national_percent': w @ p,
                       'null_national_percent': w @ n,
                       'national_model_sse': (w @ (y - p)) ** 2,
                       'national_null_sse': (w @ (y - n)) ** 2})
    annual = pd.DataFrame(annual)
    summary = []
    for (analysis, country), z in annual.groupby(['analysis', 'country'], sort=True):
        rng = np.random.default_rng(20261014)
        n = len(z)
        starts = rng.integers(0, n, size=(1999, int(np.ceil(n / 3))))
        ids = ((starts[:, :, None] + np.arange(3)) % n).reshape(1999, -1)[:, :n]
        for scale in ['field', 'national']:
            m = z[f'{scale}_model_sse'].to_numpy(); null = z[f'{scale}_null_sse'].to_numpy()
            boot = 100 * (1 - m[ids].sum(axis=1) / null[ids].sum(axis=1))
            lo, hi = np.quantile(boot, [.025, .975])
            summary.append({'analysis': analysis, 'country': country, 'scale': scale,
                            'skill_percent': 100 * (1 - m.sum() / null.sum()),
                            'ci_lower_percent': lo, 'ci_upper_percent': hi,
                            'years': n, 'start_year': int(z.analysis_year.min()),
                            'end_year': int(z.analysis_year.max()),
                            'minimum_regions': int(z.regions.min()),
                            'maximum_regions': int(z.regions.max()),
                            'baseline_mse': null.mean(), 'model_mse': m.mean()})
    return pd.DataFrame(summary), annual


def main():
    raw = pd.read_csv(INPUT / 'yield_qc.csv', usecols=['country', 'spatial_unit_id',
                         'analysis_year', 'yield_harvested_t_ha', 'harvested_area_ha'])
    archived = pd.read_csv(INPUT / 'published_field_oof.csv')
    check('raw_units_Argentina', raw[raw.country == 'Argentina'].spatial_unit_id.nunique() == 263)
    check('raw_units_Brazil', raw[raw.country == 'Brazil'].spatial_unit_id.nunique() == 1883)
    all_predictions = []
    replay_report = []
    for country in ['Argentina', 'Brazil']:
        print('Preparing', country, flush=True)
        panel = CountryPanel(country, raw, archived)
        for mode in ['published_replay', 'linear_loo_control', 'linear_nested']:
            print('Running', country, mode, flush=True)
            p = blocked(panel, mode)
            all_predictions.append(p)
            if mode == 'published_replay':
                z = p.merge(panel.archived, on=['country', 'region_key', 'analysis_year'],
                            suffixes=('_new', '_archived'), validate='one_to_one')
                diff = np.max(np.abs(z.model_prediction_percent_new - z.model_prediction_percent_archived))
                ndiff = np.max(np.abs(z.null_prediction_percent_new - z.null_prediction_percent_archived))
                check(f'published_replay_predictions_{country}', diff < 1e-7, float(diff))
                check(f'published_replay_null_{country}', ndiff < 1e-9, float(ndiff))
                check(f'published_replay_rows_{country}', len(z) == len(panel.archived), len(z))
                replay_report.append({'country': country, 'maximum_prediction_difference': diff,
                                      'maximum_null_difference': ndiff, 'rows': len(z)})
        print('Running', country, 'linear_forward', flush=True)
        p = forward(panel)
        check(f'forward_30_years_{country}', p.analysis_year.nunique() == 30)
        check(f'forward_no_future_training_{country}',
              bool((p.maximum_training_year <= p.analysis_year - 2).all()))
        # Perturb every held-out 2010 yield: past-derived transforms, training
        # records and reference weights must remain exactly unchanged.
        tr = tuple(range(1974, 2009))
        before = panel.prepare(tr, 'linear_forward', True).copy()
        saved = panel.Y[:, YEARS == 2010].copy()
        panel.Y[:, YEARS == 2010] += 10.
        panel.prepare.cache_clear()
        after = panel.prepare(tr, 'linear_forward', True).copy()
        btrain = before[before.analysis_year.isin(tr)].sort_values(['region_key','analysis_year'])
        atrain = after[after.analysis_year.isin(tr)].sort_values(['region_key','analysis_year'])
        check(f'held_out_yield_invariance_{country}',
              np.array_equal(btrain.yield_anomaly_percent.to_numpy(), atrain.yield_anomaly_percent.to_numpy()))
        check(f'held_out_weight_invariance_{country}',
              np.array_equal(before.reference_area_ha.to_numpy(), after.reference_area_ha.to_numpy()))
        panel.Y[:, YEARS == 2010] = saved
        panel.prepare.cache_clear()
        all_predictions.append(p)
    predictions = pd.concat(all_predictions, ignore_index=True)
    check('finite_paired_predictions', np.isfinite(predictions[['yield_anomaly_percent',
          'model_prediction_percent','null_prediction_percent','score_weight']].to_numpy()).all())
    check('unique_output_keys', not predictions.duplicated(['analysis','country','region_key','analysis_year']).any())
    summary, annual = summarize(predictions)
    predictions.to_csv(OUTPUT / 'held_out_predictions.csv', index=False)
    summary.to_csv(OUTPUT / 'validation_summary.csv', index=False)
    annual.to_csv(OUTPUT / 'annual_error_components.csv', index=False)
    pd.DataFrame(TUNING).to_csv(OUTPUT / 'inner_tuning.csv', index=False)
    pd.DataFrame(AUDIT).to_csv(OUTPUT / 'split_audit.csv', index=False)
    pd.DataFrame(CHECKS).to_csv(OUTPUT / 'checks.csv', index=False)
    pd.DataFrame(replay_report).to_csv(OUTPUT / 'published_replay_comparison.csv', index=False)
    environment = {'python': sys.version, 'platform': platform.platform(),
                   'numpy': np.__version__, 'pandas': pd.__version__,
                   'checks_passed': sum(x['passed'] for x in CHECKS),
                   'checks_failed': sum(not x['passed'] for x in CHECKS),
                   'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in INPUT.iterdir() if p.is_file()}}
    (OUTPUT / 'environment.json').write_text(json.dumps(environment, indent=2))
    print(summary.to_string(index=False), flush=True)
    print(json.dumps({'checks_passed': environment['checks_passed'], 'checks_failed': environment['checks_failed']}), flush=True)


if __name__ == '__main__':
    main()
