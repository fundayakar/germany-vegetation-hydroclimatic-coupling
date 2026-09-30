# Archived Google Earth Engine CSV exports

| Archive | Member names | Expected contents |
|---|---|---|
| `monthly_state_2000_2024.zip` | `Germany_monthly_state_YYYY.csv` | 25 annual files for 2000–2024; scripts retain cropland, grassland and forest records. |
| `august_5km_grid_2003_2018_2022.zip` | `Germany_monthly_grid_YYYY_m8.csv` | Three August event files with actual 5 km EPSG:3035 cell IDs (`i,j`), generated with a 500 m within-cell reduction. |
| `corine_forest_events_2003_2018_2022.zip` | `Germany_CORINEstable_foresttypes_monthly_YYYY.csv` | Three event years of April–October forest-class state summaries at a 500 m reduction scale. |

The raw exported CSV columns include `system:index` and an empty `.geo` placeholder from Earth Engine. The scripts read numeric fields and `region_id`; the empty geometries are not used to draw the maps. The grid cell coordinates are recovered from the EPSG:3035 cell index, while mapped cell NDVI values come from the CSVs.

The annual state ZIP combines the user's 19-file Drive download with six separately exported years. The values were not recalculated when repackaging; each CSV was copied intact. No source raster is included. The minimum information needed to regenerate each source export is in `gee/`.
