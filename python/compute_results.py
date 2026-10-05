"""Reproduce the descriptive tables for the 2026 TAC revision.

Run from any directory: python python/compute_results.py
Inputs are the three documented ZIP archives in data/raw.
"""
from pathlib import Path
from io import BytesIO
from zipfile import ZipFile
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
UPLOAD = ROOT / 'data' / 'raw'
OUT = ROOT / 'data' / 'derived'
OUT.mkdir(exist_ok=True)

def read_archive_csv(archive, member):
    with ZipFile(UPLOAD / archive) as z:
        return pd.read_csv(BytesIO(z.read(member)))

with ZipFile(UPLOAD / 'monthly_state_2000_2024.zip') as z:
    names = sorted(n for n in z.namelist() if n.endswith('.csv'))
    assert len(names) == 25, f'Expected 25 annual state exports, found {len(names)}'
    parts = [pd.read_csv(BytesIO(z.read(n))) for n in names]
panel = pd.concat(parts, ignore_index=True)
panel = panel[panel.lc_name.isin(['cropland', 'grassland', 'forest'])]
assert len(panel) == 7525
assert not panel.duplicated(['year', 'month', 'region_id', 'lc_name']).any()
assert sorted(panel.year.unique()) == list(range(2000, 2025))

def corr(x, y):
    return float(np.corrcoef(x, y)[0, 1])

def block_corr(df, nboot=5000, seed=1945):
    years = np.sort(df.year.unique())
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(len(years), [1/len(years)]*len(years), size=nboot)
    rows = []
    for month in [4, 5, 6, 7, 8, 9, 10]:
        dat = df[df.month.eq(month)]
        for variable in ['sm_shallow_anom', 'sm_rootzone_anom', 'vpd_anom_kpa']:
            classes = ['cropland', 'grassland', 'forest']
            observed, replicates = {}, {}
            for c in classes:
                q=dat[dat.lc_name.eq(c)]
                observed[c]=corr(q[variable],q.ndvi_anom)
                rows_year=[]
                for year in years:
                    s=q[q.year.eq(year)]
                    x=s[variable].to_numpy(); y=s.ndvi_anom.to_numpy()
                    rows_year.append([len(x),x.sum(),y.sum(),(x*x).sum(),(y*y).sum(),(x*y).sum()])
                sums=counts@np.asarray(rows_year)
                n,sx,sy,sxx,syy,sxy=sums.T
                replicates[c]=(sxy-sx*sy/n)/np.sqrt((sxx-sx*sx/n)*(syy-sy*sy/n))
            for c in classes:
                lo, hi = np.percentile(replicates[c], [2.5, 97.5])
                rows.append(dict(month=month, predictor=variable, land_cover=c, n=len(dat[dat.lc_name.eq(c)]), r=observed[c], ci_low=lo, ci_high=hi))
            diff = replicates['cropland'] - replicates['forest']
            lo, hi = np.percentile(diff, [2.5, 97.5])
            rows.append(dict(month=month, predictor=variable, land_cover='cropland_minus_forest', n=np.nan,
                             r=observed['cropland']-observed['forest'], ci_low=lo, ci_high=hi))
    return pd.DataFrame(rows)

print('Computing year-block bootstrap, 5,000 replicates per month and variable...', flush=True)
correlations = block_corr(panel)
correlations.to_csv(OUT / 'monthly_year_block_correlations.csv', index=False)

# Remove a separate linear time trend within each state, land cover and month.
cols=['ndvi_anom','sm_shallow_anom','sm_rootzone_anom','vpd_anom_kpa']
detr=panel.copy()
for key, indices in detr.groupby(['region_id','lc_name','month']).groups.items():
    q=detr.loc[indices]
    t=q.year.to_numpy(dtype=float); t=t-t.mean()
    for col in cols:
        y=q[col].to_numpy(dtype=float)
        detr.loc[indices,col]=y-(y.mean()+t*np.dot(t,y)/np.dot(t,t))
detr_corr=block_corr(detr)
detr_corr.to_csv(OUT/'monthly_detrended_year_block_correlations.csv',index=False)

events = panel[panel.year.isin([2003, 2018, 2022])].groupby(['year','month','lc_name']).agg(
    ndvi_mean=('ndvi_anom','mean'), ndvi_median=('ndvi_anom','median'),
    n_states=('ndvi_anom','size'), negative_states=('ndvi_anom',lambda x:int((x<0).sum())),
    shallow_mean=('sm_shallow_anom','mean'),root_mean=('sm_rootzone_anom','mean'),
    vpd_mean=('vpd_anom_kpa','mean')).reset_index()
events.to_csv(OUT / 'monthly_event_state_means.csv',index=False)

sub = pd.concat([read_archive_csv('corine_forest_events_2003_2018_2022.zip',
                                 f'Germany_CORINEstable_foresttypes_monthly_{y}.csv')
                 for y in [2003,2018,2022]], ignore_index=True)
assert len(sub) == 1302 and not sub.duplicated(['year','month','region_id','lc_name']).any()
subsummary = sub.groupby(['year','month','lc_name']).agg(
    ndvi_mean=('ndvi_anom','mean'), n_states=('ndvi_anom','size'),
    negative_states=('ndvi_anom',lambda x:int((x<0).sum())),
    median_fraction=('lc_fraction','median')).reset_index()
subsummary.to_csv(OUT / 'corine_event_monthly_state_means.csv',index=False)

pairs=[]
for y in [2003,2018,2022]:
    s=sub[(sub.year.eq(y)) & (sub.month.eq(8)) & (sub.lc_name.isin(['broadleaf_clc','coniferous_clc']))]
    p=s.pivot(index='region_id',columns='lc_name',values=['ndvi_anom','lc_fraction']).dropna()
    for cutoff in [0,0.01]:
        q=p[(p['lc_fraction'].min(axis=1)>=cutoff)]
        d=q['ndvi_anom']['broadleaf_clc']-q['ndvi_anom']['coniferous_clc']
        pairs.append(dict(year=y,month=8,minimum_fraction=cutoff,n_states=len(q),mean_paired_difference=d.mean(),
                          broadleaf_less_negative=int((d>0).sum())))
pd.DataFrame(pairs).to_csv(OUT/'corine_august_paired_comparison.csv',index=False)

spatial=[]
for year in [2003,2018,2022]:
    g=read_archive_csv('august_5km_grid_2003_2018_2022.zip',
                       f'Germany_monthly_grid_{year}_m8.csv')
    for lc in ['cropland','grassland','forest']:
        c=g[g.lc_name.eq(lc)].copy()
        for thresh in [0,0.2]:
            q=c[c.lc_fraction.ge(thresh)]
            w=q.lc_fraction.to_numpy()
            a=q.ndvi_anom.to_numpy()
            spatial.append(dict(year=year,class_name=lc,min_fraction=thresh,n_cells=len(q),
                sampled_fraction_sum=w.sum(),weighted_mean=np.average(a,weights=w),
                weighted_negative_share=w[a<0].sum()/w.sum()))
pd.DataFrame(spatial).to_csv(OUT/'august_grid_sensitivity.csv',index=False)

print('August correlations:')
print(correlations[correlations.month.eq(8)].to_string(index=False, float_format=lambda n:f'{n:.3f}'))
print('August CORINE pairs:')
print(pd.DataFrame(pairs).to_string(index=False,float_format=lambda n:f'{n:.4f}'))
print('Grid:')
print(pd.DataFrame(spatial).to_string(index=False,float_format=lambda n:f'{n:.3f}'))
