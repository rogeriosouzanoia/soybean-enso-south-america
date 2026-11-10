"""Behavioral checks for the numerical revision code; no synthetic paper results."""
import numpy as np
import pandas as pd
from revision_models import episode_folds, features, summarize, validate
from run_water_memory import bucket, stages, et0
from run_climate_revision import transform


def main():
    # Exact daily mass balance and storage bounds in a dry-down/recharge case.
    rain=np.r_[np.zeros(25),200.,np.zeros(25)];demand=np.full(len(rain),6.)
    s,a=bucket(rain,demand);before=np.r_[100.,s[:-1]]
    drainage=np.maximum(0,before+rain-100)
    assert np.allclose(before+rain-a-drainage,s)
    assert np.all((s>=0)&(s<=100)) and np.all(a<=demand)
    kc,middle=stages(151);assert middle.sum()==61 and np.allclose(kc[middle],1.15)
    # A sufficiently long wet period erases differences in initial storage.
    s1,_=bucket(np.full(40,10.),np.full(40,4.),initial=100.)
    s2,_=bucket(np.full(40,10.),np.full(40,4.),initial=50.)
    assert np.isclose(s1[-1],s2[-1])
    # No storage reset at sowing in a dry continuous simulation.
    s,_=bucket(np.zeros(200),np.full(200,4.));assert s[99]<1.
    # Train-only transformations cannot be changed by held-out values.
    d=pd.DataFrame({'unit':['A']*20,'year':np.arange(2000,2020),'raw':np.arange(20.)+.5})
    x=transform(d,['unit'],'raw',tuple(range(2000,2015)))
    d.loc[d.year>=2015,'raw']+=100
    z=transform(d,['unit'],'raw',tuple(range(2000,2015)))
    assert np.allclose(x[x.year<2015].value,z[z.year<2015].value)
    # Perfectly persistent lag supplies no residual innovation to the model.
    tr=pd.DataFrame({'year':np.arange(10),'oni':np.linspace(-2,2,10)})
    tr['lag']=tr.oni
    f,_=features(tr,tr,'ENSO_lag');assert np.max(np.abs(f[:,2:]))<1e-12
    # Episode allocation never separates consecutive same-phase years.
    phases={y:('La Nina' if y%5 in [0,1] else 'Neutral') for y in range(1985,2025)}
    folds=episode_folds(range(1985,2025),phases,5)
    for y in range(1986,2025):
        if phases[y]==phases[y-1]:assert folds[y]==folds[y-1]
    print('Passed: water conservation, storage bounds, stage allocation, initialization convergence, no sowing reset, training-only transformation, persistence residualization, episode integrity.')


if __name__=='__main__':main()
