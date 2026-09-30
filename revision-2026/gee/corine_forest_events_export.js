/* Monthly CORINE forest-type sensitivity for drought years in Germany.
 * Three export tasks: 2003, 2018, and 2022; April to October per task.
 * Forest pixels meet stable broad MODIS 2001-2023 mask; CORINE subtype
 * 311/312/313 agrees across its 2000/2006/2012/2018 inventories.
 */

var YEARS_TO_EXPORT = [2003, 2018, 2022];
var MONTHS = [4, 5, 6, 7, 8, 9, 10];
var START_YEAR = 2000;
var END_YEAR = 2024;
var CLASSES = ['forest_modis', 'broadleaf_clc',
               'coniferous_clc', 'mixed_clc'];
var FOLDER = 'Germany_VegClimate_revision';

var germany = ee.FeatureCollection('FAO/GAUL/2015/level0')
  .filter(ee.Filter.eq('ADM0_NAME', 'Germany')).geometry();
var states = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM0_NAME', 'Germany'))
  .map(function(f) {
    return f.set('region_id', f.get('ADM1_NAME'));
  });

// Keep the broad-class mask used for the main monthly state panel.
// The end date is exclusive, so this spans 2001-2023.
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

// CORINE records broadleaf (311), coniferous (312), and mixed (313).
// Retain only locations with the same CORINE class in all four inventories.
// Intersect with the long-term stable MODIS broad forest mask.
var clcYears = [2000, 2006, 2012, 2018];
var clcAnnual = ee.ImageCollection.fromImages(clcYears.map(function(y) {
  return ee.Image('COPERNICUS/CORINE/V20/100m/' + y)
    .select('landcover');
}));
var clcMin = clcAnnual.min();
var clcMax = clcAnnual.max();
var clcStable = clcMin.eq(clcMax);
var clcForestType = clcMin
  .updateMask(clcStable)
  .updateMask(broadLC.eq(1));

function maskForClass(name) {
  if (name === 'forest_modis') return broadLC.eq(1);
  if (name === 'broadleaf_clc') return clcForestType.eq(311);
  if (name === 'coniferous_clc') return clcForestType.eq(312);
  if (name === 'mixed_clc') return clcForestType.eq(313);
  throw new Error('Unknown forest class: ' + name);
}

var modis = ee.ImageCollection('MODIS/061/MOD13Q1')
  .filterBounds(germany)
  .filterDate('2000-01-01', '2025-01-01')
  .select('NDVI')
  .map(function(img) {
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
  return ee.Image.cat([ndvi, temp, vpd, shallow, root, precip]);
}

var YEARS = ee.List.sequence(START_YEAR, END_YEAR);
function monthlyAnomaly(y, m) {
  var observation = monthlyImage(y, m);
  var baseline = ee.ImageCollection.fromImages(YEARS.map(function(otherYear) {
    return monthlyImage(ee.Number(otherYear), m);
  })).mean();
  return observation.subtract(baseline)
    .rename(['ndvi_anom', 'temp_anom_c', 'vpd_anom_kpa',
             'sm_shallow_anom', 'sm_rootzone_anom', 'precip_anom_mm'])
    .addBands(observation.rename(['ndvi', 'temp_c', 'vpd_kpa',
                                  'sm_shallow', 'sm_rootzone', 'precip_mm']));
}

function exportClassMonth(y, m, name, monthImage) {
  var mask = maskForClass(name);
  var img = monthImage.updateMask(mask)
    .addBands(mask.unmask(0).rename('lc_fraction'));
  var result = img.reduceRegions({
    collection: states,
    reducer: ee.Reducer.mean(),
    scale: 500,
    tileScale: 16
  });
  return result.map(function(f) {
    return ee.Feature(null, {
      year: y, month: m, mode: 'state', reduction_scale_m: 500,
      region_id: f.get('region_id'), lc_name: name,
      ndvi_anom: f.get('ndvi_anom'), ndvi: f.get('ndvi'),
      temp_anom_c: f.get('temp_anom_c'), temp_c: f.get('temp_c'),
      vpd_anom_kpa: f.get('vpd_anom_kpa'), vpd_kpa: f.get('vpd_kpa'),
      sm_shallow_anom: f.get('sm_shallow_anom'),
      sm_rootzone_anom: f.get('sm_rootzone_anom'),
      sm_shallow: f.get('sm_shallow'), sm_rootzone: f.get('sm_rootzone'),
      precip_anom_mm: f.get('precip_anom_mm'),
      precip_mm: f.get('precip_mm'), lc_fraction: f.get('lc_fraction')
    });
  });
}

YEARS_TO_EXPORT.forEach(function(year) {
  var pieces = [];
  MONTHS.forEach(function(m) {
    var monthImage = monthlyAnomaly(year, m);
    CLASSES.forEach(function(name) {
      pieces.push(exportClassMonth(year, m, name, monthImage));
    });
  });
  var output = ee.FeatureCollection(pieces).flatten()
    .filter(ee.Filter.notNull(['ndvi_anom']));
  var description = 'Germany_CORINEstable_foresttypes_monthly_' + year;
  Export.table.toDrive({
    collection: output,
    description: description,
    folder: FOLDER,
    fileNamePrefix: description,
    fileFormat: 'CSV'
  });
});

Map.centerObject(germany, 6);
Map.addLayer(germany, {color: 'black'}, 'Germany AOI', false);
print('Three exports are configured for:', YEARS_TO_EXPORT);
print('Start each task in the Tasks panel.');
