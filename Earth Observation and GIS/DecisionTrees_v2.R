#==================================================================================
# R script for practical class Image Analysis
#Last modified: 03.05.2021
#==================================================================================
#Download/add the required packages/libraries
library (raster)
library (rgdal)
library(rpart)
library(rpart.plot)
library(rasterVis)
library(dismo)
library(caret)


#==================================================================================
#Define the working directory
setwd ('your working directory')
#==================================================================================
#Load the land cover data (polygon shapefile)
samples <- readOGR('Flevoland_Samples_2021.shp')

# Generate  training samples points using Flevoland_Samples dataset
pt_samples <- spsample(samples, 3000, type='regular')

# Add the land cover class to the points  
pt_samples$Class <- over(pt_samples, samples)$Class
pt_samples$Class

#writeOGR(pt_samples, dsn='your folder', layer = 'TS_Flevoland', driver="ESRI Shapefile")
#==================================================================================
#Load the Planet image and rename the spectral bands
planet_image <- stack('Flevoland_planet.tif')
names(planet_image) <- c('blue', 'green', 'red', 'NIR')
#==================================================================================
# Extracting spectral information from Planet images for the training samples generated in the previous step
sample_values <- extract(planet_image, pt_samples, df = TRUE)
sample_values <- sample_values[, -1]


# Combine the class information with extracted values
sample_data <- data.frame(classvalue = pt_samples@data$Class, sample_values)

#Split the samples into training, testing and validation  samples

Samples_partition <- createDataPartition(y = sample_data$classvalue, p = 0.70, list = FALSE)
forTraining <- sample_data[Samples_partition,]
validation <- sample_data[-Samples_partition,]


# Create training and testing datasets
Training_partition <- createDataPartition(y = forTraining$classvalue, p = 0.60, list = FALSE)
training <- forTraining[Training_partition,]
testing <- forTraining[-Training_partition,]
table(training$classvalue)
table(testing$classvalue)
table(validation$classvalue)


#==================================================================================
#Decision Tree /CART using rpart package
# Step 1: trainign classification model
DT_training <- rpart(classvalue~., data=training, method = 'class',  minsplit=yourparameter, parms = list(split = 'information'))

#Plot the trained decision tree model
rpart.plot (DT_training, type=1,  digits=3, extra=101,  fallen.leaves = TRUE)

#Investigate the trained decision tree model
#summary (DT_training)
#alternatively you can use print(DT_training)

#Apply the fitted Decision Trees model to predict the testing samples
DT_prediction <- predict(DT_training, testing, type='class')

# Assess the accuracy of the testing samples
confusionMatrix(table(DT_prediction, testing$classvalue))
confusionMatrix(table(DT_prediction, testing$classvalue))$byClass[, 1]


#investigate the complexity table. This table will be used to prune the trained DT model
printcp(DT_training)
plotcp(DT_training)
min_cp = which.min(DT_training$cptable[,"xerror"])
min_cp
#Pruning training DT model
pruned_DT <- prune(DT_training, cp= yourCPValue) 


#Predict the testing samples using the pruned decision trees and assessing the results

DT_prediction_pruned <- predict(pruned_DT, testing, type='class')

confusionMatrix(table(DT_prediction_pruned, testing$classvalue))

confusionMatrix(table(DT_prediction_pruned, testing$classvalue))$byClass[, 1]

#==================================================================================
# Step 2: Classify the entire study area
DT_prediction <- predict(planet_image, pruned_DT, type='class')
DT_prediction 
#==================================================================================

#Display the classification results
DT_prediction <- ratify(DT_prediction)
LC_raster_classes <- levels(DT_prediction)[[1]]
planetclass <- c("Arable", "BuiltUp", "Forest", "Grassland", "Water")
LC_Classes <- data.frame(classvalue1 = c(1,2,3,4,5), classnames1 = planetclass)
LC_ClassesColor <- c("#FBF65D", "#D2CDC0", "#38814E", "#D1D182", "#5475A8")
LC_raster_classes$legend <- LC_Classes$classnames
levels(DT_prediction) <- LC_raster_classes
levelplot(DT_prediction, maxpixels = 1e6,
          col.regions = LC_ClassesColor,
          scales=list(draw=FALSE),
          main = "Decision Tree classification of land cover using Planet Images")


# You can export the classification results and investigate their accuracy outside Rstudio


