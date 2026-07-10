# Germany NDVI Hydroclimatic Coupling

Analysis scripts and processed datasets for a study examining 
land-cover-specific hydroclimatic coupling of vegetation anomalies 
across Germany during the MODIS era (2000–2024).

## Repository structure
germany-ndvi-hydroclimatic-coupling/
│
├── README.md
├── LICENSE
│
├── gee/
│   └── germany_modis_v2.js          # GEE extraction script
│
├── analysis/
│   └── analysis_v2.py               # Ana Python analiz scripti
│
├── data/
│   └── Germany_VegClimate_v2_2000_2024.csv   # İşlenmiş panel veri
│
└── outputs/
    ├── 01_descriptive_by_landcover.csv
    ├── 02_correlation_by_landcover.csv
    ├── 03_interaction_models.csv
    ├── 04_moderated_mediation.csv
    ├── 05_drought_year_table.csv
    ├── 06_forest_lag_tests.csv
    └── 07_predictive_benchmark.csv

## Data sources

- MODIS NDVI: MOD13Q1 Collection 6.1 — https://doi.org/10.5067/MODIS/MOD13Q1.061
- ERA5-Land: Muñoz-Sabater et al. (2021) — https://doi.org/10.24381/cds.e2161bac
- Land cover: MCD12Q1 Collection 6.1 — https://doi.org/10.5067/MODIS/MCD12Q1.006
- ESA WorldCover v200: https://doi.org/10.5281/zenodo.7254221

## Requirements

**Python:**
pandas
numpy
scipy
statsmodels
scikit-learn

**Google Earth Engine:** GEE account required to run the extraction 
script. Processed data is provided in `data/` so re-extraction is 
not necessary to reproduce the analysis.

## How to reproduce

1. Run `analysis/analysis_v2.py` with the processed dataset in `data/`
2. Outputs will be saved to `outputs/`
3. GEE script in `gee/` can be used to re-extract raw data if needed

## Citation

If you use this code or data, please cite the archived version:

> Yakar, F. (2025). germany-ndvi-hydroclimatic-coupling [Software]. 
