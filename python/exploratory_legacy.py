"""Exploratory prior-summer to next-spring forest NDVI associations.

Run: python python/exploratory_legacy.py
The year-clustered models are descriptive and do not establish legacy effects.
"""
from pathlib import Path
from io import BytesIO
from zipfile import ZipFile
import numpy as np
import pandas as pd
from scipy.stats import t as t_dist

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data'/'raw'/'monthly_state_2000_2024.zip'
OUT = ROOT/'data'/'derived'
OUT.mkdir(exist_ok=True)
with ZipFile(DATA) as z:
    parts = [pd.read_csv(BytesIO(z.read(n))) for n in sorted(z.namelist()) if n.endswith('.csv')]
panel = pd.concat(parts, ignore_index=True).query("lc_name == 'forest'")
assert len(panel) == 25*7*15

cols = ['ndvi_anom','sm_rootzone_anom','vpd_anom_kpa']
def season(months, prefix):
    a = panel[panel.month.isin(months)].groupby(['region_id','year'],as_index=False)[cols].mean()
    return a.rename(columns={c:prefix+c for c in cols})

spring = season([4,5,6], 'spring_')
prior = season([6,7,8], 'prior_')
prior['year'] += 1
data = spring.merge(prior,on=['region_id','year'],validate='one_to_one')
assert len(data)==24*15

def fit(outcome, target, controls):
    predictors=[target]+controls
    q=data[['region_id','year',outcome]+predictors].dropna().copy()
    # Match the originally reported sample-standardized coefficients exactly.
    y=q[outcome].to_numpy()/q[outcome].std()
    X=np.column_stack([np.ones(len(q)),
        *(q[c].to_numpy()/q[c].std() for c in predictors),
        (q.year.to_numpy()-q.year.mean())/q.year.std(),
        pd.get_dummies(q.region_id,drop_first=True).to_numpy().astype(float)])
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    residual=y-X@beta
    inverse=np.linalg.pinv(X.T@X)
    years=q.year.to_numpy();meat=np.zeros((X.shape[1],X.shape[1]))
    for year in np.unique(years):
        score=X[years==year].T@residual[years==year]
        meat+=np.outer(score,score)
    G=q.year.nunique();N=len(q);K=X.shape[1]
    variance=inverse@meat@inverse*G/(G-1)*(N-1)/(N-K)
    se=float(np.sqrt(max(0,variance[1,1])))
    estimate=float(beta[1])
    return dict(outcome=outcome,previous_summer_predictor=target,
                controls=', '.join(controls),n_state_years=N,n_years=G,
                standardized_beta=estimate,year_cluster_se=se,
                two_sided_p=float(2*t_dist.sf(abs(estimate/se),G-1)))

controls=['spring_sm_rootzone_anom','spring_vpd_anom_kpa','prior_ndvi_anom']
results=pd.DataFrame([fit('spring_ndvi_anom',target,controls)
                      for target in ['prior_sm_rootzone_anom','prior_vpd_anom_kpa']])
results.to_csv(OUT/'exploratory_next_spring_forest.csv',index=False)
print(results.to_string(index=False,float_format=lambda x:f'{x:.3f}'))
