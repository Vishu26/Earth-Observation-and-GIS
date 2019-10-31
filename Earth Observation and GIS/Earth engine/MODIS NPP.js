var dataset = ee.ImageCollection('MODIS/006/MOD17A3H')
                  .filterDate('2017-10-01','2018-10-28');

var meanIm = dataset.select('Npp').mean();
var fapar = ee.Image(meanIm);
var geometry = ee.Geometry.Rectangle([72.623387, 23.184179, 72.633213,
23.192804]);
var nppVis = {
  min: 0.0,
  max: 19000.0,
  palette: ['bbe029', '0a9501', '074b03'],
};
//Map.setCenter(geometry, 13)
Map.addLayer(fapar.clip(geometry), nppVis)
Export.image.toDrive({
  image: fapar,
  description: 'imageToDriveExample1',
  scale: 10,
  region: geometry
});