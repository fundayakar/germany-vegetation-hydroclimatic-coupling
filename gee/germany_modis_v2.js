/***************************************
 GERMANY MODIS VEGETATION-CLIMATE v2
 Period: 2000-2024
 New: VPD, precipitation, root-zone SM,
      antecedent SM, stable pixel mask,
      year-by-year export (memory safe)
***************************************/

var START_YEAR = 2000;
var END_YEAR = 2024;
var EXPORT_FOLDER = 'GEE_Germany_VegClimate_v2';

// ---- EXPORT ONE YEAR AT A TIME ----
// Change this to export different years.
// Run once per year, or use the batch
// loop at the bottom (may hit memory).
var EXPORT_YEAR = 2000; // <-- change this manually per run

var years = ee.List.sequence(START_YEAR, END_YEAR);

// ---- GEOMETRY ----
var germany = ee.FeatureCollection('FAO/GAUL/2015/level0')
  .filter(ee.Filter.eq('ADM0_NAME', 'Germany'))
  .geometry();

var states = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM0_NAME', 'Germany'));


/***************************************
 STABLE PIXEL MASK (MODIS MCD12Q1)
 Keeps only pixels that stayed in the
 same LC class throughout 2001-2023.
 Uses IGBP classification:
   Forest    = classes 1-5
   Grassland = class 10
   Cropland  = classes 12, 14
***************************************/
var mcd12 = ee.ImageCollection('MODIS/061/MCD12Q1')
  .filterBounds(germany)
  .filterDate('2001-01-01', '2024-01-01')
  .select('LC_Type1');

// Reclassify each annual image to 3 classes
// 1 = forest, 2 = grassland, 3 = cropland, 0 = other
function reclassMCD12(img) {
  var forest    = img.lte(5).and(img.gte(1));
  var grassland = img.eq(10);
  var cropland  = img.eq(12).or(img.eq(14));
  return forest.multiply(1)
    .add(grassland.multiply(2))
    .add(cropland.multiply(3))
    .rename('lc_reclass')
    .copyProperties(img, ['system:time_start']);
}

var mcd12Reclass = mcd12.map(reclassMCD12);

// Mode across all years = dominant class
var lcMode = mcd12Reclass.reduce(ee.Reducer.mode()).rename('lc_mode');

// Stability: std dev across years = 0 means never changed
var lcStd = mcd12Reclass.reduce(ee.Reducer.stdDev()).rename('lc_std');

// Stable mask: same class every year AND not "other"
var stableMask = lcStd.eq(0).and(lcMode.neq(0));
var stableLC   = lcMode.updateMask(stableMask).clip(germany);

// Also keep ESA WorldCover for comparison/sensitivity
var esa = ee.Image('ESA/WorldCover/v200/2021').select('Map').clip(germany);
var esaReclass = esa.remap([10, 30, 40], [1, 2, 3]).rename('lc_esa');
esaReclass = esaReclass.updateMask(esaReclass.neq(0));

// PRIMARY LC source for analysis: stableLC (MCD12Q1)
// SENSITIVITY LC source: esaReclass (ESA WorldCover)
// Change lcSource below to switch:
var lcSource = stableLC; // or: esaReclass


/***************************************
 MODIS NDVI (MOD13Q1, 250m, 16-day)
 Available from 2000-02-18
***************************************/
var modis = ee.ImageCollection('MODIS/061/MOD13Q1')
  .filterBounds(germany)
  .filterDate(
    ee.Date.fromYMD(START_YEAR, 1, 1),
    ee.Date.fromYMD(END_YEAR + 1, 1, 1)
  )
  .map(function(img) {
    return img.select('NDVI')
      .multiply(0.0001)
      .rename('NDVI')
      .copyProperties(img, ['system:time_start']);
  });


/***************************************
 ERA5-LAND (monthly)
 New variables added:
  - total_precipitation_sum  (precip)
  - dewpoint_temperature_2m  (for VPD)
  - volumetric_soil_water_layer_3 (root-zone)
  - volumetric_soil_water_layer_4 (deep)
***************************************/
var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
  .filterBounds(germany)
  .filterDate(
    ee.Date.fromYMD(START_YEAR - 1, 1, 1), // -1 year for antecedent
    ee.Date.fromYMD(END_YEAR + 1, 1, 1)
  )
  .map(function(img) {
    var t2m = img.select('temperature_2m')
      .subtract(273.15)
      .rename('t2m_c');

    var td = img.select('dewpoint_temperature_2m')
      .subtract(273.15)
      .rename('td_c');

    // VPD = saturation vapour pressure - actual vapour pressure
    // es = 0.6108 * exp(17.27 * T / (T + 237.3))  [kPa]
    // ea = 0.6108 * exp(17.27 * Td / (Td + 237.3)) [kPa]
    // VPD = es - ea
    var es = img.select('temperature_2m').subtract(273.15)
      .expression(
        '0.6108 * exp(17.27 * T / (T + 237.3))',
        { T: img.select('temperature_2m').subtract(273.15) }
      ).rename('es_kpa');

    var ea = img.select('dewpoint_temperature_2m').subtract(273.15)
      .expression(
        '0.6108 * exp(17.27 * Td / (Td + 237.3))',
        { Td: img.select('dewpoint_temperature_2m').subtract(273.15) }
      ).rename('ea_kpa');

    var vpd = es.subtract(ea).rename('vpd_kpa');

    // Soil moisture layers
    var sm1 = img.select('volumetric_soil_water_layer_1').rename('sm_layer1'); // 0-7cm
    var sm2 = img.select('volumetric_soil_water_layer_2').rename('sm_layer2'); // 7-28cm
    var sm3 = img.select('volumetric_soil_water_layer_3').rename('sm_layer3'); // 28-100cm
    var sm4 = img.select('volumetric_soil_water_layer_4').rename('sm_layer4'); // 100-289cm

    // Shallow SM (surface, for cropland/grassland): layers 1+2
    var smShallow = sm1.add(sm2).divide(2).rename('sm_shallow');

    // Root-zone SM (for forests): layers 3+4 depth-weighted
    // layer3 = 72cm depth, layer4 = 189cm depth → weights ~0.38 and 0.62
    var smRootzone = sm3.multiply(0.38).add(sm4.multiply(0.62)).rename('sm_rootzone');

    // Profile SM (all layers combined): weighted by layer thickness
    // layer1=7cm, layer2=21cm, layer3=72cm, layer4=189cm → total=289cm
    var smProfile = sm1.multiply(7).add(sm2.multiply(21))
      .add(sm3.multiply(72)).add(sm4.multiply(189))
      .divide(289).rename('sm_profile');

    // Precipitation (monthly sum, m → mm)
    var precip = img.select('total_precipitation_sum')
      .multiply(1000)
      .rename('precip_mm');

    return t2m
      .addBands(td)
      .addBands(vpd)
      .addBands(sm1)
      .addBands(sm2)
      .addBands(sm3)
      .addBands(sm4)
      .addBands(smShallow)
      .addBands(smRootzone)
      .addBands(smProfile)
      .addBands(precip)
      .copyProperties(img, ['system:time_start']);
  });


/***************************************
 COMPUTE LONG-TERM MEANS (for anomalies)
 Baseline: 2000-2024 for climate
           2000-2024 for NDVI
 Using full period mean as baseline.
 (Can switch to 2001-2020 if preferred)
***************************************/

// NDVI: growing season mean per year
var annualNDVI = ee.ImageCollection.fromImages(
  years.map(function(y) {
    y = ee.Number(y);
    var growingStart = ee.Date.fromYMD(y, 4, 1);
    var growingEnd   = ee.Date.fromYMD(y, 11, 1);
    var summerStart  = ee.Date.fromYMD(y, 6, 1);
    var summerEnd    = ee.Date.fromYMD(y, 9, 1);

    var growing = modis.filterDate(growingStart, growingEnd).mean()
      .rename('growing_mean_ndvi');
    var summer  = modis.filterDate(summerStart, summerEnd).mean()
      .rename('summer_mean_ndvi');

    return growing.addBands(summer)
      .clip(germany)
      .set('year', y)
      .set('system:time_start', ee.Date.fromYMD(y, 1, 1).millis());
  })
);

var ltGrowing = annualNDVI.select('growing_mean_ndvi').mean();
var ltSummer  = annualNDVI.select('summer_mean_ndvi').mean();

// Climate: summer means per year
// Summer = June-August (JJA)
// Spring = March-May (MAM) for antecedent SM
var annualClimate = ee.ImageCollection.fromImages(
  years.map(function(y) {
    y = ee.Number(y);

    var summerStart  = ee.Date.fromYMD(y, 6, 1);
    var summerEnd    = ee.Date.fromYMD(y, 9, 1);
    var springStart  = ee.Date.fromYMD(y, 3, 1);
    var springEnd    = ee.Date.fromYMD(y, 6, 1);
    // Antecedent: previous year summer precip
    var prevSumStart = ee.Date.fromYMD(y.subtract(1), 6, 1);
    var prevSumEnd   = ee.Date.fromYMD(y.subtract(1), 9, 1);

    var summer = era5.filterDate(summerStart, summerEnd);
    var spring = era5.filterDate(springStart, springEnd);
    var prevSum = era5.filterDate(prevSumStart, prevSumEnd);

    return summer.select('t2m_c').mean().rename('summer_temp')
      .addBands(summer.select('vpd_kpa').mean().rename('summer_vpd'))
      .addBands(summer.select('sm_shallow').mean().rename('summer_sm_shallow'))
      .addBands(summer.select('sm_rootzone').mean().rename('summer_sm_rootzone'))
      .addBands(summer.select('sm_profile').mean().rename('summer_sm_profile'))
      .addBands(summer.select('precip_mm').sum().rename('summer_precip'))
      // Spring SM: proxy for pre-season water availability
      .addBands(spring.select('sm_shallow').mean().rename('spring_sm_shallow'))
      .addBands(spring.select('sm_rootzone').mean().rename('spring_sm_rootzone'))
      // Antecedent summer precip (previous year)
      .addBands(prevSum.select('precip_mm').sum().rename('prev_summer_precip'))
      .clip(germany)
      .set('year', y)
      .set('system:time_start', ee.Date.fromYMD(y, 1, 1).millis());
  })
);

// Long-term means for anomaly calculation
var ltTemp          = annualClimate.select('summer_temp').mean();
var ltVPD           = annualClimate.select('summer_vpd').mean();
var ltSmShallow     = annualClimate.select('summer_sm_shallow').mean();
var ltSmRootzone    = annualClimate.select('summer_sm_rootzone').mean();
var ltSmProfile     = annualClimate.select('summer_sm_profile').mean();
var ltPrecip        = annualClimate.select('summer_precip').mean();
var ltSpringShallow = annualClimate.select('spring_sm_shallow').mean();
var ltSpringRoot    = annualClimate.select('spring_sm_rootzone').mean();
var ltPrevPrecip    = annualClimate.select('prev_summer_precip').mean();


/***************************************
 BUILD ANOMALIES FOR A SINGLE YEAR
 Called once per export task
***************************************/
function buildAnomalyImage(y) {
  y = ee.Number(y);

  var ndvi = ee.Image(
    annualNDVI.filter(ee.Filter.eq('year', y)).first()
  );
  var clim = ee.Image(
    annualClimate.filter(ee.Filter.eq('year', y)).first()
  );

  // NDVI anomalies
  var ndviAnom = ndvi.select('growing_mean_ndvi')
    .subtract(ltGrowing).rename('anom_growing_ndvi');
  var ndviSumAnom = ndvi.select('summer_mean_ndvi')
    .subtract(ltSummer).rename('anom_summer_ndvi');

  // Climate anomalies
  var tempAnom     = clim.select('summer_temp')
    .subtract(ltTemp).rename('anom_temp');
  var vpdAnom      = clim.select('summer_vpd')
    .subtract(ltVPD).rename('anom_vpd');
  var smShallowAnom = clim.select('summer_sm_shallow')
    .subtract(ltSmShallow).rename('anom_sm_shallow');
  var smRootzoneAnom = clim.select('summer_sm_rootzone')
    .subtract(ltSmRootzone).rename('anom_sm_rootzone');
  var smProfileAnom = clim.select('summer_sm_profile')
    .subtract(ltSmProfile).rename('anom_sm_profile');
  var precipAnom   = clim.select('summer_precip')
    .subtract(ltPrecip).rename('anom_precip');
  var springShallowAnom = clim.select('spring_sm_shallow')
    .subtract(ltSpringShallow).rename('anom_spring_sm_shallow');
  var springRootAnom = clim.select('spring_sm_rootzone')
    .subtract(ltSpringRoot).rename('anom_spring_sm_rootzone');
  var prevPrecipAnom = clim.select('prev_summer_precip')
    .subtract(ltPrevPrecip).rename('anom_prev_precip');

  return ndvi
    .addBands(clim)
    .addBands(ndviAnom)
    .addBands(ndviSumAnom)
    .addBands(tempAnom)
    .addBands(vpdAnom)
    .addBands(smShallowAnom)
    .addBands(smRootzoneAnom)
    .addBands(smProfileAnom)
    .addBands(precipAnom)
    .addBands(springShallowAnom)
    .addBands(springRootAnom)
    .addBands(prevPrecipAnom)
    .set('year', y);
}


/***************************************
 ZONAL STATISTICS FOR ONE LC CLASS
***************************************/
function statsForLC(img, lcVal, lcName) {
  var masked = img.updateMask(lcSource.eq(lcVal));

  var stats = masked.reduceRegions({
    collection: states,
    reducer: ee.Reducer.mean(),
    scale: 5000,
    tileScale: 16  // increased from 8 to reduce memory
  });

  return stats.map(function(f) {
    return ee.Feature(null, {
      year:                    img.get('year'),
      adm1_name:               f.get('ADM1_NAME'),
      lc_class:                lcVal,
      lc_name:                 lcName,
      // Raw NDVI
      growing_mean_ndvi:       f.get('growing_mean_ndvi'),
      summer_mean_ndvi:        f.get('summer_mean_ndvi'),
      // NDVI anomalies
      anom_growing_ndvi:       f.get('anom_growing_ndvi'),
      anom_summer_ndvi:        f.get('anom_summer_ndvi'),
      // Raw climate
      summer_temp:             f.get('summer_temp'),
      summer_vpd:              f.get('summer_vpd'),
      summer_sm_shallow:       f.get('summer_sm_shallow'),
      summer_sm_rootzone:      f.get('summer_sm_rootzone'),
      summer_sm_profile:       f.get('summer_sm_profile'),
      summer_precip:           f.get('summer_precip'),
      spring_sm_shallow:       f.get('spring_sm_shallow'),
      spring_sm_rootzone:      f.get('spring_sm_rootzone'),
      prev_summer_precip:      f.get('prev_summer_precip'),
      // Climate anomalies
      anom_temp:               f.get('anom_temp'),
      anom_vpd:                f.get('anom_vpd'),
      anom_sm_shallow:         f.get('anom_sm_shallow'),
      anom_sm_rootzone:        f.get('anom_sm_rootzone'),
      anom_sm_profile:         f.get('anom_sm_profile'),
      anom_precip:             f.get('anom_precip'),
      anom_spring_sm_shallow:  f.get('anom_spring_sm_shallow'),
      anom_spring_sm_rootzone: f.get('anom_spring_sm_rootzone'),
      anom_prev_precip:        f.get('anom_prev_precip')
    });
  });
}


/***************************************
 EXPORT: SINGLE YEAR MODE (recommended)

 Instructions:
 1. Set EXPORT_YEAR at the top of script
 2. Click Run → Tasks tab → Run the task
 3. Change EXPORT_YEAR, repeat for each year
 4. Merge all CSVs in Python (see note below)

 Python merge (after all exports):
   import pandas as pd, glob
   files = glob.glob('Germany_VegClimate_v2_*.csv')
   df = pd.concat([pd.read_csv(f) for f in files])
   df.to_csv('Germany_VegClimate_v2_2000_2024.csv', index=False)
***************************************/

var img = buildAnomalyImage(EXPORT_YEAR);

var zonalTable = ee.FeatureCollection([
  statsForLC(img, 1, 'forest'),
  statsForLC(img, 2, 'grassland'),
  statsForLC(img, 3, 'cropland')
]).flatten();

Export.table.toDrive({
  collection: zonalTable,
  description: 'Germany_VegClimate_v2_' + EXPORT_YEAR,
  folder: EXPORT_FOLDER,
  fileNamePrefix: 'Germany_VegClimate_v2_' + EXPORT_YEAR,
  fileFormat: 'CSV'
});


/***************************************
 OPTIONAL: BATCH EXPORT ALL YEARS
 Uncomment this block instead of the
 single-year export above if your GEE
 account has enough quota.
 Risk: may hit memory on some years.
 Safer alternative: use single-year mode.

var yearList = ee.List.sequence(START_YEAR, END_YEAR);
yearList.evaluate(function(yList) {
  yList.forEach(function(y) {
    var img = buildAnomalyImage(y);
    var zonalTable = ee.FeatureCollection([
      statsForLC(img, 1, 'forest'),
      statsForLC(img, 2, 'grassland'),
      statsForLC(img, 3, 'cropland')
    ]).flatten();
    Export.table.toDrive({
      collection: zonalTable,
      description: 'Germany_VegClimate_v2_' + y,
      folder: EXPORT_FOLDER,
      fileNamePrefix: 'Germany_VegClimate_v2_' + y,
      fileFormat: 'CSV'
    });
  });
});
***************************************/


/***************************************
 QUICK DIAGNOSTIC (optional)
 Print first year to Console to check
 values before running full export.
 Uncomment to use.

var testImg = buildAnomalyImage(2018); // drought year
print('2018 anomaly image bands:', testImg.bandNames());

var testStats = statsForLC(testImg, 3, 'cropland')
  .limit(3);
print('Cropland sample 2018:', testStats);
***************************************/
