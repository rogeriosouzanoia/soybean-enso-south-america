from pathlib import Path
import csv
from decimal import Decimal
from itertools import groupby
from collections import Counter

ROOT=Path(__file__).resolve().parent
def read(name):
    with (ROOT/name).open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
s=sorted(read("ENSO_ONI_seasonal.csv"),key=lambda r:r["center_date"])
annual=read("ENSO_indices_by_soybean_harvest_year_1974_2024.csv")
def candidate(r):
    x=Decimal(r["oni_c"])
    return "El Nino" if x>=Decimal(".50") else "La Nina" if x<=Decimal("-.50") else "Neutral"
bydate={}
for run,(sign,records) in enumerate(groupby(s,key=candidate),1):
    rows=list(records)
    phase=sign if sign!="Neutral" and len(rows)>=5 else "Neutral"
    peak=max(abs(Decimal(r["oni_c"])) for r in rows) if phase!="Neutral" else None
    strength="Neutral" if peak is None else "Weak" if peak<1 else "Moderate" if peak<Decimal("1.50") else "Strong" if peak<2 else "Very strong"
    peaks=[r["season"]+" "+r["year"] for r in rows if peak is not None and abs(Decimal(r["oni_c"]))==peak]
    for r in rows:
        assert r["oni_event_phase"]==phase
        assert r["oni_event_strength"]==strength
        assert int(r["oni_run_length"])==len(rows)
        if peak is not None:assert Decimal(r["oni_event_peak_absolute_c"])==peak
        bydate[r["center_date"]]=(r,run,phase,strength,peak,rows[0]["center_date"],rows[-1]["center_date"],peaks)
out=[]
for a in annual:
    y=int(a["harvest_year"]);r,run,phase,strength,peak,start,end,peaks=bydate[f"{y}-01-01"]
    assert a["oni_harvest_year_phase"]==phase
    assert a["oni_harvest_year_strength"]==strength
    assert Decimal(a["oni_djf_c"])==Decimal(r["oni_c"])
    out.append(dict(harvest_year=y,oni_djf_c=r["oni_c"],phase=phase,episode_strength=strength,display_strength="Strong" if strength=="Very strong" else strength,qualifying_run_id=run if peak is not None else "",episode_start_center_date=start if peak is not None else "",episode_end_center_date=end if peak is not None else "",episode_peak_absolute_c=str(peak) if peak is not None else "",peak_seasons="; ".join(peaks)))
assert len(out)==51 and len({r["harvest_year"] for r in out})==51
with (ROOT/"Annual_ENSO_episode_classification.csv").open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
counts=Counter(r["phase"] for r in out)
assert counts=={"El Nino":17,"La Nina":19,"Neutral":15}
print("Verified all seasonal assignments and 51 harvest years:",dict(counts))
