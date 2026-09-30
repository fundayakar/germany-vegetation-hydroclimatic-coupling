# Germany vegetation and hydroclimatic coupling

Code and processed datasets for Funda Yakar's study of land-cover-specific vegetation anomalies across Germany during 2000–2024.

## Analysis versions

- **2026 major revision:** [`revision-2026/`](revision-2026/) contains the monthly state panel, actual 5 km August grid checks, CORINE forest-class diagnostics, complete extraction and analysis scripts, seven numerical tables and five figures. Its [README](revision-2026/README.md) gives the run order, input provenance and methodological limits. This work is under manuscript revision; figures and claims may still change before journal resubmission.
- **Earlier annual analysis:** [`gee/germany_modis_v2.js`](gee/germany_modis_v2.js), [`python/analysis_v2_final.py`](python/analysis_v2_final.py), [`Germany_VegClimate_v2_2000_2024.csv`](Germany_VegClimate_v2_2000_2024.csv) and [`outputs/`](outputs/) are retained for provenance. Their methods and conclusions should not be substituted for the monthly revision.

The original state extraction used an Earth Engine reduction scale of 5000 m rather than a common 5 km grid. The new grid is constructed explicitly in the revision folder. The revision also documents the April–October monthly windows, retrospective land-cover masks and limitations of MODIS NDVI and ERA5-Land.

## Primary data

- MODIS MOD13Q1 Collection 6.1 NDVI: https://doi.org/10.5067/MODIS/MOD13Q1.061
- MODIS MCD12Q1 Collection 6.1 land cover: https://doi.org/10.5067/MODIS/MCD12Q1.061
- ERA5-Land: https://doi.org/10.5194/essd-13-4349-2021
- CORINE Land Cover: https://land.copernicus.eu/en/products/corine-land-cover
- FAO GAUL 2015 administrative polygons, accessed through Google Earth Engine

The earlier GEE script also contains a WorldCover comparison. The primary class masks for both the annual and monthly analyses use MCD12Q1.

## Software licence

The existing [`LICENSE`](LICENSE) applies to repository software. Source datasets remain subject to their providers' terms. A citable manuscript version and any archival DOI should be added only after the revised manuscript and bibliography are finalized.
