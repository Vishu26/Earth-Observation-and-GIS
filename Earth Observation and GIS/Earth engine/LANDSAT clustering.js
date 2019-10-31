// Load a pre-computed Landsat composite for input.
var input = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filter(ee.Filter.date('2018-01-01', '2018-05-01'));

var mean = input.select(['B2', 'B3', 'B4', 'B5']).mean()
var im = ee.Image(mean)
// Define a region in which to generate a sample of the input.
var region = ee.Geometry.Rectangle(29.7, 30, 32.5, 31.7);

// Display the sample region.
Map.setCenter(73.6834, 20.7500,15);
Map.addLayer(ee.Image().paint(region, 0, 2), {}, 'region');

// Make the training dataset.
var training = im.sample({
  region: region,
  scale: 30,
  numPixels: 5000
});

// Instantiate the clusterer and train it.
var clusterer = ee.Clusterer.wekaKMeans(7).train(training);

// Cluster the input using the trained clusterer.
var result = im.cluster(clusterer);

// Display the clusters with random colors.
Map.addLayer(result.randomVisualizer(), {}, 'clusters');