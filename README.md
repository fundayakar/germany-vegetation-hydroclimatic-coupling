# Germany vegetation and hydroclimate revision analyses

Code, processed Earth Engine exports and figures for the revision of **“Land-Cover-Specific Hydroclimatic Coupling of Vegetation Anomalies across Germany during the MODIS Era”** (*Theoretical and Applied Climatology*). This repository contains the 2026 monthly and spatial sensitivity analysis. The manuscript and reviewer correspondence are maintained separately.

## Contents

| Path | Contents |
|---|---|
| `gee/monthly_landcover_export.js` | Full Google Earth Engine script for monthly broad-class state exports and actual 5 km August grid exports. Configure `YEAR`, `MODE` and `GRID_MONTH` at the top before each task. |
| `gee/corine_forest_events_export.js` | Full Earth Engine script for April–October CORINE forest-class exports in 2003, 2018 and 2022. It defines three tasks. |
| `data/raw/*.zip` | The 25 annual monthly state CSVs, three August grid CSVs and three CORINE forest-class CSVs as originally exported. See `data/raw/README.md`. |
| `python/compute_results.py` | Monthly Pearson correlations with year-block intervals, event means, grid sensitivity and paired forest-class contrasts. |
| `python/exploratory_legacy.py` | Exploratory prior-summer to next-spring forest regressions with year-clustered standard errors. |
| `python/make_figures.py` | Generates the five manuscript figures from the archived inputs and derived tables. |
| `data/derived/*.csv` | Recomputed numerical outputs, including unrounded values. |
| `figures/*.png` | Five figures numbered as in the revised manuscript: study area (1), monthly coupling (2), event trajectories (3), August grid maps (4), and forest classes (5). |

## Reproduce the analysis

Python 3.11 or newer is recommended. The workflow was checked using Python 3.12, numpy 2.3.5, pandas 2.2.3, scipy 1.17.0 and matplotlib 3.10.8. From **this directory**:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python python/compute_results.py
python python/exploratory_legacy.py
python python/make_figures.py
```

The scripts read directly from the three ZIP archives. No manual decompression or Earth Engine login is needed to reproduce the tables and figures. The scripts overwrite the corresponding files under `data/derived` and `figures`. To regenerate the primary exports from satellite and reanalysis imagery, use the full GEE scripts and their Google Earth Engine dataset IDs; exports must be started in the Tasks panel. The 2003 CORINE job was resource-intensive and may time out if quotas differ.

## What the analyses compare

The main panel contains 25 years (2000–2024), seven calendar months (April–October) and 43 available state/class combinations per month: 15 cropland, 13 grassland and 15 forest, totalling 7525 observations. The original Earth Engine state extraction used `reduceRegions(scale: 5000)`; this is a **state summary at a 5000 m reduction scale**, not an intermediate 5 km grid. Calendar-month anomalies use each month's 2000–2024 pixel baseline. The year-block bootstrap resamples the 25 years jointly across states and classes (5000 draws, fixed seed 1945). Detrended correlations first remove a separate linear year trend within each state, class and calendar month.

The actual grid analysis uses 5 km equal-area EPSG:3035 cells in three August events, with image reduction at 500 m within each cell. The `lc_fraction` weight represents the fraction of *sampled stable-class locations* in a cell. It is not a validated fraction of all forest land. Grid maps and weighted cell means should not be numerically compared with equal-state means as if only the aggregation unit changed; their reduction scales and weights also differ. ERA5-Land hydroclimate remains much coarser than the MODIS samples and cells.

For the forest composition check, CORINE 311/312/313 classifications must agree across four inventories (2000, 2006, 2012 and 2018) and overlap the stable MODIS broad forest mask. MCD12Q1 broad-class stability is evaluated over 2001–2023 because the GEE filter end date `2024-01-01` is exclusive. CORINE is a mapped land-cover product (100 m raster, approximately 25 ha mapping unit), so its classes are not pure species stands. The mixed class has little sampled coverage, and a ≥1% sampled state-share check is included for the broad-leaved versus coniferous comparison.

The MOD13Q1 NDVI extraction matches the earlier analysis and does not apply an additional quality-assurance screen. The retrospective stable-class mask, NDVI saturation and artifacts, year dependence, event selection and correlated hydroclimatic predictors limit inference. The legacy regressions are descriptive rather than causal evidence of delayed physiological effects.

## Provenance and reuse

Primary inputs are MOD13Q1 and MCD12Q1 Collection 6.1 (NASA LP DAAC), ERA5-Land Monthly Aggregated (Copernicus Climate Data Store) and CORINE Land Cover (Copernicus Land Monitoring Service), processed using FAO GAUL 2015 Germany and state boundaries in Earth Engine. Archive ZIP checksums are in `data/raw/SHA256SUMS.txt`. These source products retain their providers' respective terms; no separate licence is asserted for them here.

The data and code availability statement in the manuscript should cite this repository and, if an archived release is deposited, its verified DOI. No Zenodo record for this revised analysis is asserted here.
