"""
Germany NDVI Hydroclimatic Coupling Analysis
Land-cover-specific hydroclimatic coupling of vegetation anomalies
across Germany during the MODIS era (2000-2024)

Outputs (saved to OUTPUT_DIR):
  01_descriptive_by_landcover.csv
  02_correlation_by_landcover.csv
  03_interaction_models.csv
  03b_standardized_coefficients.csv
  03c_marginal_vpd_slopes.csv
  03d_clustered_se_robustness.csv
  04_moderated_mediation.csv
  05_drought_year_table.csv
  06_forest_lag_tests.csv
  07_predictive_benchmark.csv

Requirements:
  pip install pandas numpy scipy statsmodels scikit-learn
"""

import pandas as pd
import numpy as np
from scipy import stats as scipy_stats
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, GroupKFold
import warnings
import os

warnings.filterwarnings('ignore')

# ── CONFIG ───────────────────────────────────────────────────────────────────
DATA_PATH   = 'Germany_VegClimate_v2_2000_2024.csv'  # update path if needed
OUTPUT_DIR  = 'outputs'
TARGET      = 'anom_growing_ndvi'
REGION_FE   = 'adm1_name'
LC_COL      = 'lc_name'
YEAR_COL    = 'year'
DROUGHT_YRS = [2003, 2018, 2019, 2022]
N_BOOT      = 1000
RANDOM_SEED = 42

os.makedirs(OUTPUT_DIR, exist_ok=True)
np.random.seed(RANDOM_SEED)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────
df = pd.read_csv(DATA_PATH)
df = df.dropna(subset=[TARGET])
df['lc_factor'] = pd.Categorical(
    df[LC_COL], categories=['cropland', 'grassland', 'forest'])
print(f"Loaded: {df.shape[0]} rows, {df.shape[1]} cols")
print(f"Years: {sorted(df[YEAR_COL].unique())}")
print(f"LC:    {df[LC_COL].value_counts().to_dict()}\n")

PRED_COLS = [
    'anom_sm_shallow', 'anom_sm_rootzone', 'anom_sm_profile',
    'anom_temp', 'anom_vpd',
    'anom_precip', 'anom_prev_precip',
    'anom_spring_sm_shallow', 'anom_spring_sm_rootzone'
]

M3_FORMULA = (
    f'{TARGET} ~ anom_sm_shallow + anom_vpd + anom_precip + anom_temp '
    f'+ anom_spring_sm_shallow '
    f'+ C({LC_COL}, Treatment("cropland")) '
    f'+ anom_sm_shallow:C({LC_COL}, Treatment("cropland")) '
    f'+ anom_vpd:C({LC_COL}, Treatment("cropland")) '
    f'+ C({REGION_FE})'
)


# ════════════════════════════════════════════════════════════════════════════
# 01  DESCRIPTIVE STATISTICS BY LAND COVER
# ════════════════════════════════════════════════════════════════════════════
print("=" * 60)
print("01  DESCRIPTIVE BY LAND COVER")
print("=" * 60)

anom_cols = [c for c in df.columns if c.startswith('anom_')]
rows = []
for lc in ['cropland', 'grassland', 'forest']:
    sub = df[df[LC_COL] == lc]
    for col in anom_cols:
        rows.append({
            'lc': lc, 'variable': col,
            'n':      sub[col].notna().sum(),
            'mean':   sub[col].mean(),
            'std':    sub[col].std(),
            'min':    sub[col].min(),
            'p25':    sub[col].quantile(0.25),
            'median': sub[col].median(),
            'p75':    sub[col].quantile(0.75),
            'max':    sub[col].max()
        })

desc_df = pd.DataFrame(rows).round(5)
desc_df.to_csv(f'{OUTPUT_DIR}/01_descriptive_by_landcover.csv', index=False)
print(desc_df[desc_df['variable'] == TARGET][
    ['lc', 'mean', 'std', 'min', 'max']].to_string(index=False))
print()


# ════════════════════════════════════════════════════════════════════════════
# 02  CORRELATIONS BY LAND COVER + VIF
# ════════════════════════════════════════════════════════════════════════════
print("=" * 60)
print("02  CORRELATIONS BY LAND COVER")
print("=" * 60)

corr_rows = []
for lc in ['cropland', 'grassland', 'forest']:
    sub = df[df[LC_COL] == lc].dropna(subset=PRED_COLS + [TARGET])
    for col in PRED_COLS:
        r, p = scipy_stats.pearsonr(sub[col], sub[TARGET])
        sig = ('***' if p < 0.001 else '**' if p < 0.01
               else '*' if p < 0.05 else 'ns')
        corr_rows.append({'lc': lc, 'predictor': col,
                          'r': round(r, 4), 'p': round(p, 4),
                          'sig': sig, 'n': len(sub)})

corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(f'{OUTPUT_DIR}/02_correlation_by_landcover.csv', index=False)

for lc in ['cropland', 'grassland', 'forest']:
    print(f"\n{lc.upper()}:")
    print(corr_df[corr_df['lc'] == lc][
        ['predictor', 'r', 'p', 'sig']].to_string(index=False))

# VIF
print("\nVIF (cropland subset, core predictors):")
vif_cols = ['anom_sm_shallow', 'anom_temp', 'anom_vpd',
            'anom_precip', 'anom_spring_sm_shallow']
sub_vif = df[df[LC_COL] == 'cropland'][vif_cols].dropna()
X_vif = sm.add_constant(sub_vif)
vif_vals = [variance_inflation_factor(X_vif.values, i)
            for i in range(X_vif.shape[1])]
vif_df = pd.DataFrame({'variable': X_vif.columns, 'VIF': vif_vals}).round(2)
print(vif_df.to_string(index=False))
print()


# ════════════════════════════════════════════════════════════════════════════
# 03  INTERACTION MODELS
# ════════════════════════════════════════════════════════════════════════════
print("=" * 60)
print("03  INTERACTION MODELS")
print("=" * 60)

FORMULAS = {
    'M1_SM_only':
        f'{TARGET} ~ anom_sm_shallow + C({LC_COL}, Treatment("cropland")) '
        f'+ anom_sm_shallow:C({LC_COL}, Treatment("cropland")) '
        f'+ C({REGION_FE})',
    'M2_atm_demand':
        f'{TARGET} ~ anom_vpd + anom_temp + anom_precip '
        f'+ C({LC_COL}, Treatment("cropland")) '
        f'+ anom_vpd:C({LC_COL}, Treatment("cropland")) '
        f'+ C({REGION_FE})',
    'M3_combined': M3_FORMULA,
    'M4_combined_rootzone':
        f'{TARGET} ~ anom_sm_shallow + anom_sm_rootzone + anom_vpd '
        f'+ anom_precip + anom_temp + anom_spring_sm_shallow '
        f'+ C({LC_COL}, Treatment("cropland")) '
        f'+ anom_sm_shallow:C({LC_COL}, Treatment("cropland")) '
        f'+ anom_sm_rootzone:C({LC_COL}, Treatment("cropland")) '
        f'+ anom_vpd:C({LC_COL}, Treatment("cropland")) '
        f'+ C({REGION_FE})',
    'M3b_yearFE_robustness':
        f'{TARGET} ~ anom_sm_shallow + anom_vpd + anom_precip + anom_temp '
        f'+ C({LC_COL}, Treatment("cropland")) '
        f'+ anom_sm_shallow:C({LC_COL}, Treatment("cropland")) '
        f'+ C({REGION_FE}) + C({YEAR_COL})',
}

model_results = []
fitted_models = {}

for name, formula in FORMULAS.items():
    try:
        mod = smf.ols(formula, data=df).fit(cov_type='HC3')
        fitted_models[name] = mod
        fe_terms = [t for t in mod.params.index
                    if 'adm1' in t.lower() or 'C(year' in t]
        core_params = mod.params.drop(fe_terms)
        core_pvals  = mod.pvalues.drop(fe_terms)
        core_ci     = mod.conf_int().drop(fe_terms)
        core_bse    = mod.bse.drop(fe_terms)

        for param in core_params.index:
            sig = ('***' if core_pvals[param] < 0.001 else
                   '**'  if core_pvals[param] < 0.01  else
                   '*'   if core_pvals[param] < 0.05  else 'ns')
            model_results.append({
                'model':   name,
                'term':    param,
                'coef':    round(core_params[param], 5),
                'se':      round(core_bse[param], 5),
                'ci_lo':   round(core_ci.loc[param, 0], 5),
                'ci_hi':   round(core_ci.loc[param, 1], 5),
                'p':       round(core_pvals[param], 4),
                'sig':     sig,
                'n':       int(mod.nobs),
                'R2':      round(mod.rsquared, 4),
                'adj_R2':  round(mod.rsquared_adj, 4)
            })
        print(f"{name}: n={int(mod.nobs)}, R²={mod.rsquared:.3f}, "
              f"adj.R²={mod.rsquared_adj:.3f}")
    except Exception as e:
        print(f"{name}: ERROR - {e}")

model_df = pd.DataFrame(model_results)
model_df.to_csv(f'{OUTPUT_DIR}/03_interaction_models.csv', index=False)
print()

# ── 03b  STANDARDIZED COEFFICIENTS (M3) ──────────────────────────────────
print("03b  STANDARDIZED COEFFICIENTS (M3)")
pred_std_cols = ['anom_sm_shallow', 'anom_vpd', 'anom_precip',
                 'anom_temp', 'anom_spring_sm_shallow']
df_std = df.copy()
for col in pred_std_cols + [TARGET]:
    df_std[col] = (df[col] - df[col].mean()) / df[col].std()

mod_std = smf.ols(M3_FORMULA, data=df_std).fit(cov_type='HC3')
param_names_std = list(mod_std.params.index)
fe_terms_std = [t for t in param_names_std if 'adm1' in t.lower()]
core_std = [(t, mod_std.params[t], mod_std.bse[t],
             mod_std.conf_int().loc[t, 0], mod_std.conf_int().loc[t, 1],
             mod_std.pvalues[t])
            for t in param_names_std if t not in fe_terms_std]

std_rows = []
for t, b, se, lo, hi, p in core_std:
    sig = ('***' if p < 0.001 else '**' if p < 0.01
           else '*' if p < 0.05 else 'ns')
    std_rows.append({'term': t, 'beta_std': round(b, 4),
                     'se': round(se, 4), 'ci_lo': round(lo, 4),
                     'ci_hi': round(hi, 4), 'p': round(p, 4), 'sig': sig,
                     'n': int(mod_std.nobs),
                     'R2': round(mod_std.rsquared, 4)})
    print(f"  {t[:50]:<50} β={b:+.3f} SE={se:.3f} p={p:.4f} {sig}")

std_df = pd.DataFrame(std_rows)
std_df.to_csv(f'{OUTPUT_DIR}/03b_standardized_coefficients.csv', index=False)
print()

# ── 03c  MARGINAL VPD SLOPES (delta method) ───────────────────────────────
print("03c  MARGINAL VPD SLOPES")
mod_m3 = fitted_models['M3_combined']
param_names = list(mod_m3.params.index)

b_vpd_crop  = mod_m3.params['anom_vpd']
se_vpd_crop = mod_m3.bse['anom_vpd']
p_vpd_crop  = mod_m3.pvalues['anom_vpd']
ci_vpd_crop = mod_m3.conf_int().loc['anom_vpd']

cov = mod_m3.cov_params()

def marginal_slope(base_term, int_term):
    b  = mod_m3.params[base_term] + mod_m3.params[int_term]
    v  = (cov.loc[base_term, base_term]
          + cov.loc[int_term, int_term]
          + 2 * cov.loc[base_term, int_term])
    se = np.sqrt(v)
    t  = b / se
    p  = 2 * scipy_stats.t.sf(abs(t), df=mod_m3.df_resid)
    ci_lo = b - 1.96 * se
    ci_hi = b + 1.96 * se
    return b, se, p, ci_lo, ci_hi

b_f, se_f, p_f, lo_f, hi_f = marginal_slope(
    'anom_vpd',
    'anom_vpd:C(lc_name, Treatment("cropland"))[T.forest]')
b_g, se_g, p_g, lo_g, hi_g = marginal_slope(
    'anom_vpd',
    'anom_vpd:C(lc_name, Treatment("cropland"))[T.grassland]')

slope_rows = [
    {'lc': 'cropland',  'slope': round(b_vpd_crop, 4),
     'se': round(se_vpd_crop, 4), 'p': round(p_vpd_crop, 4),
     'ci_lo': round(ci_vpd_crop[0], 4), 'ci_hi': round(ci_vpd_crop[1], 4)},
    {'lc': 'forest',    'slope': round(b_f, 4), 'se': round(se_f, 4),
     'p': round(p_f, 4), 'ci_lo': round(lo_f, 4), 'ci_hi': round(hi_f, 4)},
    {'lc': 'grassland', 'slope': round(b_g, 4), 'se': round(se_g, 4),
     'p': round(p_g, 4), 'ci_lo': round(lo_g, 4), 'ci_hi': round(hi_g, 4)},
]
slope_df = pd.DataFrame(slope_rows)
slope_df.to_csv(f'{OUTPUT_DIR}/03c_marginal_vpd_slopes.csv', index=False)
for r in slope_rows:
    print(f"  {r['lc']:<12} slope={r['slope']:+.4f} SE={r['se']:.4f} "
          f"p={r['p']:.4f} 95%CI=[{r['ci_lo']:.4f}, {r['ci_hi']:.4f}]")
print()

# ── 03d  CLUSTERED SE ROBUSTNESS ──────────────────────────────────────────
print("03d  CLUSTERED SE ROBUSTNESS")
mod_ols = smf.ols(M3_FORMULA, data=df).fit()
param_names_cl = list(mod_ols.params.index)

hc3_res   = mod_ols.get_robustcov_results(cov_type='HC3')
state_res = mod_ols.get_robustcov_results(
    cov_type='cluster', groups=df[REGION_FE].values)

hc3_p   = pd.Series(hc3_res.pvalues,   index=param_names_cl)
state_p = pd.Series(state_res.pvalues, index=param_names_cl)

df['state_year'] = df[REGION_FE].astype(str) + '_' + df[YEAR_COL].astype(str)
V_state = np.array(state_res.cov_params())
V_year  = np.array(mod_ols.get_robustcov_results(
    cov_type='cluster', groups=df[YEAR_COL].values).cov_params())
V_sy    = np.array(mod_ols.get_robustcov_results(
    cov_type='cluster', groups=df['state_year'].values).cov_params())
V_tw    = V_state + V_year - V_sy
diag_tw = np.diag(V_tw)
diag_tw = np.where(diag_tw < 0, 0, diag_tw)
se_tw   = np.sqrt(diag_tw)

cl_terms = [
    'anom_sm_shallow', 'anom_vpd', 'anom_precip', 'anom_temp',
    'anom_vpd:C(lc_name, Treatment("cropland"))[T.forest]',
    'anom_vpd:C(lc_name, Treatment("cropland"))[T.grassland]'
]

cl_rows = []
print(f"  {'Term':<50} {'HC3_p':>7} {'State_p':>9} {'2way_p':>8}")
print("  " + "-"*76)
for t in cl_terms:
    if t not in param_names_cl:
        continue
    idx   = param_names_cl.index(t)
    b     = mod_ols.params[t]
    p_h   = hc3_p[t]
    p_s   = state_p[t]
    se_t  = se_tw[idx]
    p_tw  = (2 * scipy_stats.t.sf(abs(b / se_t), df=mod_ols.df_resid)
             if se_t > 0 else np.nan)
    label = t[:50]
    print(f"  {label:<50} {p_h:>7.4f} {p_s:>9.4f} "
          f"{p_tw:>8.4f}" if not np.isnan(p_tw) else
          f"  {label:<50} {p_h:>7.4f} {p_s:>9.4f} {'NA':>8}")
    cl_rows.append({'term': t, 'coef': round(b, 5),
                    'p_hc3': round(p_h, 4), 'p_state_cl': round(p_s, 4),
                    'p_twoway_cl': round(p_tw, 4) if not np.isnan(p_tw) else None})

cl_df = pd.DataFrame(cl_rows)
cl_df.to_csv(f'{OUTPUT_DIR}/03d_clustered_se_robustness.csv', index=False)
print()


# ════════════════════════════════════════════════════════════════════════════
# 04  MODERATED MEDIATION
# ════════════════════════════════════════════════════════════════════════════
print("=" * 60)
print("04  MODERATED MEDIATION")
print("=" * 60)

med_rows = []

for lc in ['cropland', 'grassland', 'forest']:
    sub = df[df[LC_COL] == lc].dropna(
        subset=['anom_temp', 'anom_vpd', 'anom_sm_shallow',
                TARGET, REGION_FE])

    # a paths
    a1_mod = smf.ols(
        f'anom_sm_shallow ~ anom_temp + C({REGION_FE})',
        data=sub).fit(cov_type='HC3')
    a1 = a1_mod.params['anom_temp']

    a2_mod = smf.ols(
        f'anom_sm_shallow ~ anom_vpd + C({REGION_FE})',
        data=sub).fit(cov_type='HC3')
    a2 = a2_mod.params['anom_vpd']

    # b path (temp mediation)
    b_mod = smf.ols(
        f'{TARGET} ~ anom_sm_shallow + anom_temp + C({REGION_FE})',
        data=sub).fit(cov_type='HC3')
    b      = b_mod.params['anom_sm_shallow']
    c_prime = b_mod.params['anom_temp']

    # b path (VPD mediation)
    b_vpd_mod = smf.ols(
        f'{TARGET} ~ anom_sm_shallow + anom_vpd + C({REGION_FE})',
        data=sub).fit(cov_type='HC3')
    b_vpd = b_vpd_mod.params['anom_sm_shallow']

    # Bootstrap
    indirect_temp_boot = []
    indirect_vpd_boot  = []
    idx_arr = np.arange(len(sub))

    for _ in range(N_BOOT):
        boot = sub.iloc[np.random.choice(idx_arr, size=len(idx_arr),
                                          replace=True)]
        try:
            a1_b = smf.ols(
                f'anom_sm_shallow ~ anom_temp + C({REGION_FE})',
                data=boot).fit().params['anom_temp']
            a2_b = smf.ols(
                f'anom_sm_shallow ~ anom_vpd + C({REGION_FE})',
                data=boot).fit().params['anom_vpd']
            b_b = smf.ols(
                f'{TARGET} ~ anom_sm_shallow + anom_temp + C({REGION_FE})',
                data=boot).fit().params['anom_sm_shallow']
            b_vpd_b = smf.ols(
                f'{TARGET} ~ anom_sm_shallow + anom_vpd + C({REGION_FE})',
                data=boot).fit().params['anom_sm_shallow']
            indirect_temp_boot.append(a1_b * b_b)
            indirect_vpd_boot.append(a2_b * b_vpd_b)
        except Exception:
            pass

    ind_temp = a1 * b
    ind_vpd  = a2 * b_vpd
    ci_temp  = (np.percentile(indirect_temp_boot, 2.5),
                np.percentile(indirect_temp_boot, 97.5))
    ci_vpd   = (np.percentile(indirect_vpd_boot, 2.5),
                np.percentile(indirect_vpd_boot, 97.5))
    sig_temp = 'sig' if (ci_temp[0] > 0 or ci_temp[1] < 0) else 'ns'
    sig_vpd  = 'sig' if (ci_vpd[0]  > 0 or ci_vpd[1]  < 0) else 'ns'

    print(f"\n{lc.upper()} (n={len(sub)}):")
    print(f"  a1 (temp→SM):    {a1:.5f}")
    print(f"  a2 (VPD→SM):     {a2:.5f}")
    print(f"  b  (SM→NDVI, temp model): {b:.5f}")
    print(f"  b  (SM→NDVI, VPD model):  {b_vpd:.5f}")
    print(f"  Indirect (temp): {ind_temp:.6f} 95%CI "
          f"[{ci_temp[0]:.6f}, {ci_temp[1]:.6f}] {sig_temp}")
    print(f"  Indirect (VPD):  {ind_vpd:.6f} 95%CI "
          f"[{ci_vpd[0]:.6f}, {ci_vpd[1]:.6f}] {sig_vpd}")

    med_rows.append({
        'lc': lc, 'n': len(sub),
        'a1_temp_to_SM':   round(a1, 5),
        'a2_vpd_to_SM':    round(a2, 5),
        'b_SM_to_NDVI_temp_model': round(b, 5),
        'b_SM_to_NDVI_vpd_model':  round(b_vpd, 5),
        'c_prime_temp':    round(c_prime, 5),
        'indirect_temp':   round(ind_temp, 6),
        'indirect_temp_ci_lo': round(ci_temp[0], 6),
        'indirect_temp_ci_hi': round(ci_temp[1], 6),
        'indirect_temp_sig': sig_temp,
        'indirect_vpd':    round(ind_vpd, 6),
        'indirect_vpd_ci_lo': round(ci_vpd[0], 6),
        'indirect_vpd_ci_hi': round(ci_vpd[1], 6),
        'indirect_vpd_sig': sig_vpd,
        'b_R2_temp': round(b_mod.rsquared, 4),
        'b_R2_vpd':  round(b_vpd_mod.rsquared, 4)
    })

med_df = pd.DataFrame(med_rows)
med_df.to_csv(f'{OUTPUT_DIR}/04_moderated_mediation.csv', index=False)
print()


# ════════════════════════════════════════════════════════════════════════════
# 05  DROUGHT-YEAR EVENT TABLE
# ════════════════════════════════════════════════════════════════════════════
print("=" * 60)
print("05  DROUGHT-YEAR EVENT TABLE")
print("=" * 60)

drought_rows = []
for y in DROUGHT_YRS:
    for lc in ['cropland', 'grassland', 'forest']:
        sub = df[(df[YEAR_COL] == y) & (df[LC_COL] == lc)]
        drought_rows.append({
            'year': y, 'lc': lc,
            'anom_growing_ndvi':  round(sub['anom_growing_ndvi'].mean(), 4),
            'anom_temp':          round(sub['anom_temp'].mean(), 4),
            'anom_vpd':           round(sub['anom_vpd'].mean(), 4),
            'anom_sm_shallow':    round(sub['anom_sm_shallow'].mean(), 4),
            'anom_sm_rootzone':   round(sub['anom_sm_rootzone'].mean(), 4),
            'anom_precip':        round(sub['anom_precip'].mean(), 4),
        })

drought_df = pd.DataFrame(drought_rows)
drought_df.to_csv(f'{OUTPUT_DIR}/05_drought_year_table.csv', index=False)
pivot = drought_df.pivot_table(
    index='year', columns='lc', values='anom_growing_ndvi').round(4)
print("NDVI anomaly by drought year and land cover:")
print(pivot.to_string())
print()


# ════════════════════════════════════════════════════════════════════════════
# 06  FOREST LAG TESTS
# ════════════════════════════════════════════════════════════════════════════
print("=" * 60)
print("06  FOREST LAG TESTS")
print("=" * 60)

forest_df = df[df[LC_COL] == 'forest'].copy()
forest_df = forest_df.sort_values([REGION_FE, YEAR_COL])
forest_df['ndvi_t1'] = forest_df.groupby(REGION_FE)[TARGET].shift(-1)
forest_df['ndvi_t2'] = forest_df.groupby(REGION_FE)[TARGET].shift(-2)

lag_rows = []
for lag_col, lag_name in [
    (TARGET,      'same-year (t)'),
    ('ndvi_t1',   'lag t+1'),
    ('ndvi_t2',   'lag t+2')
]:
    for pred in ['anom_sm_shallow', 'anom_sm_rootzone',
                 'anom_temp', 'anom_vpd', 'anom_precip']:
        sub = forest_df.dropna(subset=[lag_col, pred])
        if len(sub) < 10:
            continue
        r, p = scipy_stats.pearsonr(sub[pred], sub[lag_col])
        sig = ('***' if p < 0.001 else '**' if p < 0.01
               else '*' if p < 0.05 else 'ns')
        lag_rows.append({'response': lag_name, 'predictor': pred,
                         'r': round(r, 4), 'p': round(p, 4),
                         'sig': sig, 'n': len(sub)})

lag_df = pd.DataFrame(lag_rows)
lag_df.to_csv(f'{OUTPUT_DIR}/06_forest_lag_tests.csv', index=False)
print(lag_df.to_string(index=False))

print("\nForest NDVI anomaly around key drought events:")
for dy in [2003, 2018]:
    print(f"\n  Drought year: {dy}")
    for yr in [dy - 1, dy, dy + 1, dy + 2]:
        val = forest_df[forest_df[YEAR_COL] == yr][TARGET].mean()
        print(f"    {yr}: {val:.4f}")
print()


# ════════════════════════════════════════════════════════════════════════════
# 07  PREDICTIVE BENCHMARK (Random Forest, year-grouped CV)
# ════════════════════════════════════════════════════════════════════════════
print("=" * 60)
print("07  PREDICTIVE BENCHMARK")
print("=" * 60)

FEATURE_COLS = [
    'anom_sm_shallow', 'anom_sm_rootzone', 'anom_sm_profile',
    'anom_temp', 'anom_vpd',
    'anom_precip', 'anom_prev_precip',
    'anom_spring_sm_shallow', 'anom_spring_sm_rootzone',
    'lc_class'
]

bench_rows = []
for lc_subset, label in [
    (None,       'all_LC'),
    ('cropland', 'cropland'),
    ('grassland','grassland'),
    ('forest',   'forest')
]:
    sub = df if lc_subset is None else df[df[LC_COL] == lc_subset]
    sub = sub.dropna(subset=FEATURE_COLS + [TARGET])
    X      = sub[FEATURE_COLS].values
    y      = sub[TARGET].values
    groups = sub[YEAR_COL].values

    rf = RandomForestRegressor(
        n_estimators=300, max_features='sqrt',
        min_samples_leaf=5, random_state=RANDOM_SEED, n_jobs=-1)

    n_splits = min(5, len(np.unique(groups)))
    gkf = GroupKFold(n_splits=n_splits)

    cv_r2   = cross_val_score(rf, X, y, cv=gkf, groups=groups, scoring='r2')
    cv_rmse = np.sqrt(-cross_val_score(
        rf, X, y, cv=gkf, groups=groups,
        scoring='neg_mean_squared_error'))

    rf.fit(X, y)
    importances = pd.Series(
        rf.feature_importances_, index=FEATURE_COLS).sort_values(ascending=False)

    bench_rows.append({
        'subset':          label,
        'n':               len(sub),
        'cv_R2_mean':      round(cv_r2.mean(), 3),
        'cv_R2_std':       round(cv_r2.std(), 3),
        'cv_RMSE_mean':    round(cv_rmse.mean(), 5),
        'cv_RMSE_std':     round(cv_rmse.std(), 5),
        'top_predictor_1': importances.index[0],
        'top_predictor_2': importances.index[1],
        'top_predictor_3': importances.index[2],
    })
    print(f"\n{label} (n={len(sub)}):")
    print(f"  CV R² = {cv_r2.mean():.3f} ± {cv_r2.std():.3f}")
    print(f"  Top predictors: {list(importances.index[:4])}")

bench_df = pd.DataFrame(bench_rows)
bench_df.to_csv(f'{OUTPUT_DIR}/07_predictive_benchmark.csv', index=False)

print()
print("=" * 60)
print("ALL DONE. Outputs saved to:", OUTPUT_DIR)
print("=" * 60)
