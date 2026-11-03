// PARAMETERS
var COLL_ID  = 'COPERNICUS/S3/OLCI';                          // Sentinel-3 optical (same platform as the thermal)
var SITES_ID = 'projects/dryve-512810/assets/icos_coords'; // table with site, lat, lon
var START    = '2016-04-01';
var END      = '2026-11-01';
var OUT_NAME = 's3_overpasses';

// ------------------------------------------------------------------
  var sites = ee.FeatureCollection(SITES_ID);

// for each station: one row per image that covers it
var overpasses = sites.map(function (s) {
  var coll = ee.ImageCollection(COLL_ID)
  .filterBounds(s.geometry())
  .filterDate(START, END);
  return coll.map(function (img) {
    return ee.Feature(null, {
      site:     s.get('site'),
      time_utc: img.date().format('YYYY-MM-dd HH:mm:ss'),    // UTC
      image_id: img.get('system:index')                       // contains S3A or S3B
    });
  });
}).flatten();

print('first rows', overpasses.limit(5));

Export.table.toDrive({
  collection:  overpasses,
  description: OUT_NAME,
  fileFormat:  'CSV',
  selectors:   ['site', 'time_utc', 'image_id']
});