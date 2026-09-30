/*
 * Germany vegetation revision: calendar-month NDVI and hydroclimate anomalies
 *
 * Run in the Google Earth Engine Code Editor. Export one year per task.
 * This is a NEW extraction, not a claim that the old CSV contains grid cells.
 * The code has been syntax checked locally; server execution and quotas must
 * be checked in Earth Engine before the resulting data are used in the paper.
 *
 * MODE = 'state' reproduces the original 5-km sampling scale for a matched
 * temporal-sensitivity test. MODE = 'grid' makes actual 5-km polygons, then
 * averages the original images at 500 m within each polygon. The grid does
 * not create new spatial resolution in ERA5-Land (~9 km).
 */

var YEAR = 2018;                   // change and rerun for each year
var MODE = 'state';                // 'state' or 'grid'
var FOLDER = 'Germany_VegClimate_revision';
var START_YEAR = 2000;
var END_YEAR = 2024;
var GRID_MONTH = 8;               // in grid mode, export ONE month per task
var MONTHS = (MODE === 'grid') ? [GRID_MONTH] : [4, 5, 6, 7, 8, 9, 10];
var GRID_SIZE_M = 5000;
var SAMPLE_SCALE_M = (MODE === 'grid') ? 500 : 5000;
var CLASSES = ['cropland', 'grassland', 'forest',
               'evergreen_needleleaf', 'deciduous_broadleaf', 'mixed_forest'];

if (YEAR < START_YEAR || YEAR > END_YEAR) {
  throw new Error('YEAR must be 2000 through 2024');
}
if (MODE !== 'state' && MODE !== 'grid') {
  throw new Error("MODE must be 'state' or 'grid'");
}

var germany = ee.FeatureCollection('FAO/GAUL/2015/level0')
  .filter(ee.Filter.eq('ADM0_NAME', 'Germany')).geometry();
var states = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM0_NAME', 'Germany'))
  .map(function (f) {
    return f.set('region_id', f.get('ADM1_NAME'));
  });

// Preserve the original 2001-2023 broad-class mask for a matched comparison.
// In the original script, '2024-01-01' is an exclusive filter end date.
var lcAnnual = ee.ImageCollection('MODIS/061/MCD12Q1')
  .filterBounds(germany)
  .filterDate('2001-01-01', '2024-01-01')
  .select('LC_Type1');

function broadClass(img) {
  var c = img.select('LC_Type1');
  return c.gte(1).and(c.lte(5)).multiply(1)
    .add(c.eq(10).multiply(2))
    .add(c.eq(12).or(c.eq(14)).multiply(3))
    .rename('broad');
}

var broadAnnual = lcAnnual.map(broadClass);
var broadStable = broadAnnual.reduce(ee.Reducer.stdDev()).eq(0);
var broadMode = broadAnnual.reduce(ee.Reducer.mode());
var broadLC = broadMode.updateMask(broadStable.and(broadMode.neq(0)));

// Subtype comparison is a stricter sample: exact IGBP class never changes
// during 2001-2023. Do not call these remote-sensing classes pure stands.
var rawMin = lcAnnual.min();
var rawMax = lcAnnual.max();
var stableSubtype = rawMin.eq(rawMax);
var subtypeLC = rawMin.updateMask(stableSubtype).updateMask(broadLC.eq(1));

function maskForClass(name) {
  if (name === 'forest') return broadLC.eq(1);
  if (name === 'grassland') return broadLC.eq(2);
  if (name === 'cropland') return broadLC.eq(3);
  if (name === 'evergreen_needleleaf') return subtypeLC.eq(1);
  if (name === 'deciduous_broadleaf') return subtypeLC.eq(4);
  if (name === 'mixed_forest') return subtypeLC.eq(5);
  throw new Error('Unknown land cover: ' + name);
}

var modis = ee.ImageCollection('MODIS/061/MOD13Q1')
  .filterBounds(germany)
  .filterDate('2000-01-01', '2025-01-01')
  .select('NDVI')
  .map(function (img) {
    // Match the original pipeline; a QA-screened run can be added separately.
    return img.multiply(0.0001).rename('ndvi')
      .copyProperties(img, ['system:time_start']);
  });

var era = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
  .filterBounds(germany)
  .filterDate('2000-01-01', '2025-01-01');

function monthlyImage(y, m) {
  var start = ee.Date.fromYMD(y, m, 1);
  var end = start.advance(1, 'month');
  var ndvi = modis.filterDate(start, end).mean().rename('ndvi');
  var climate = ee.Image(era.filterDate(start, end).first());
  var temp = climate.select('temperature_2m').subtract(273.15)
    .rename('temp_c');
  var dew = climate.select('dewpoint_temperature_2m').subtract(273.15);
  // Same numerical saturation-vapour-pressure approximation as the old code.
  var es = temp.expression('0.6108 * exp(17.27 * t / (t + 237.3))',
                           {t: temp});
  var ea = dew.expression('0.6108 * exp(17.27 * t / (t + 237.3))',
                          {t: dew});
  var vpd = es.subtract(ea).rename('vpd_kpa');
  var shallow = climate.select('volumetric_soil_water_layer_1')
    .add(climate.select('volumetric_soil_water_layer_2'))
    .divide(2).rename('sm_shallow');
  var root = climate.select('volumetric_soil_water_layer_3')
    .multiply(0.38)
    .add(climate.select('volumetric_soil_water_layer_4').multiply(0.62))
    .rename('sm_rootzone');
  var precip = climate.select('total_precipitation_sum')
    .multiply(1000).rename('precip_mm');
  return ee.Image.cat([ndvi, temp, vpd, shallow, root, precip])
    .set('year', y).set('month', m);
}

// Baseline is calculated separately for each calendar month, at each pixel.
// Therefore August 2018 is compared with August 2000-2024, not with April.
var YEARS = ee.List.sequence(START_YEAR, END_YEAR);
function monthlyAnomaly(y, m) {
  var observation = monthlyImage(y, m);
  var baseline = ee.ImageCollection.fromImages(YEARS.map(function (otherYear) {
    return monthlyImage(ee.Number(otherYear), m);
  })).mean();
  var names = ['ndvi_anom', 'temp_anom_c', 'vpd_anom_kpa',
               'sm_shallow_anom', 'sm_rootzone_anom', 'precip_anom_mm'];
  return observation.subtract(baseline).rename(names)
    .addBands(observation.rename(['ndvi', 'temp_c', 'vpd_kpa',
                                  'sm_shallow', 'sm_rootzone', 'precip_mm']));
}

// A real equal-area 5-km grid is created only in grid mode.
var regions;
if (MODE === 'grid') {
  regions = germany.coveringGrid(ee.Projection('EPSG:3035'), GRID_SIZE_M)
    .map(function (f) {
      return f.set('region_id', f.id());
    });
} else {
  regions = states;
}

function exportClassMonth(y, m, name, monthImage) {
  var mask = maskForClass(name);
  var img = monthImage.updateMask(mask);
  // Fraction of sampled land-cover pixels in each polygon. This is a
  // coverage diagnostic, not an exact 250-m land-cover area measurement.
  img = img.addBands(mask.unmask(0).rename('lc_fraction'));
  var reduceOptions = {
    collection: regions,
    reducer: ee.Reducer.mean(),
    scale: SAMPLE_SCALE_M,
    tileScale: 16
  };
  // Original extraction did not set a CRS. Keep that convention for the
  // state-scale temporal check; use the equal-area CRS for actual grid cells.
  if (MODE === 'grid') reduceOptions.crs = 'EPSG:3035';
  var result = img.reduceRegions(reduceOptions);
  return result.map(function (f) {
    return ee.Feature(null, {
      year: y, month: m, mode: MODE, reduction_scale_m: SAMPLE_SCALE_M,
      region_id: f.get('region_id'), lc_name: name,
      ndvi_anom: f.get('ndvi_anom'), ndvi: f.get('ndvi'),
      temp_anom_c: f.get('temp_anom_c'), temp_c: f.get('temp_c'),
      vpd_anom_kpa: f.get('vpd_anom_kpa'), vpd_kpa: f.get('vpd_kpa'),
      sm_shallow_anom: f.get('sm_shallow_anom'),
      sm_rootzone_anom: f.get('sm_rootzone_anom'),
      sm_shallow: f.get('sm_shallow'),
      sm_rootzone: f.get('sm_rootzone'),
      precip_anom_mm: f.get('precip_anom_mm'),
      precip_mm: f.get('precip_mm'),
      lc_fraction: f.get('lc_fraction')
    });
  });
}

var pieces = [];
MONTHS.forEach(function (m) {
  var monthImage = monthlyAnomaly(YEAR, m);
  CLASSES.forEach(function (name) {
    pieces.push(exportClassMonth(YEAR, m, name, monthImage));
  });
});
var output = ee.FeatureCollection(pieces).flatten()
  .filter(ee.Filter.notNull(['ndvi_anom']));

var suffix = (MODE === 'grid') ? '_m' + GRID_MONTH : '';
var description = 'Germany_monthly_' + MODE + '_' + YEAR + suffix;
Export.table.toDrive({
  collection: output,
  description: description,
  folder: FOLDER,
  fileNamePrefix: description,
  fileFormat: 'CSV'
});

Map.centerObject(germany, 6);
Map.addLayer(germany, {color: 'black'}, 'Germany AOI', false);
print('Mode:', MODE, 'Year:', YEAR, 'Months:', MONTHS);
print('One export task is configured; start it in the Tasks panel.');
