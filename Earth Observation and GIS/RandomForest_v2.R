# If not installed, you need to install all relevant packages first
#install.packages("raster", "randomForest", "sp", "rgdal", "ggplot2")

#load libraries
library(raster)
library(randomForest)
library(sp)
library(rgdal)
library(ggplot2)
library(caret)
set.seed(42)

# Setting-up the working directory 
setwd("Define your working directory")
###### Import images
S2_image  = "S2_June22_2020.tif"
S2        = stack(S2_image)
names(S2) = c('blue', 'green', 'red', 'RE_5','RE_6', 'RE_7', 'NIR', 'NIR_2', 'SWIR_1', 'SWIR_2')

#==================================================================================
#Calculate vegetation indices
#==================================================================================
NDVI_B8B5  = (S2$NIR-S2$RE_5)/(S2$NIR+S2$RE_5)
NDVI_B8B6  = (S2$NIR-S2$RE_6)/(S2$NIR+S2$RE_6)
NDVI_B8B7  = (S2$NIR-S2$RE_7)/(S2$NIR+S2$RE_7)
NDVI_B8aB5 = (S2$NIR_2-S2$RE_5)/(S2$NIR_2+S2$RE_5)
NDVI_B8aB6 = (S2$NIR_2-S2$RE_6)/(S2$NIR_2+S2$RE_6)
NDVI_B8aB7 = (S2$NIR_2-S2$RE_7)/(S2$NIR_2+S2$RE_7)
PSRI       = (S2$NIR-S2$green)/S2$RE_5
Clre       = (S2$RE_7/S2$RE_5)/1
NDre1      =  (S2$RE_6-S2$RE_5)/(S2$RE_6+S2$RE_5)
NDre2      =  (S2$RE_7-S2$RE_5)/(S2$RE_7+S2$RE_5)
NDWI_1     =  (S2$green-S2$NIR)/(S2$green+S2$NIR)
NDWI_2     =  (S2$green-S2$SWIR_1)/(S2$green+S2$SWIR_1)  
NDWI_3     =  (S2$green-S2$SWIR_2)/(S2$green+S2$SWIR_2) 

S2_blue    = S2$blue/10000
S2_green   = S2$green/10000
S2_red     = S2$red/10000
S2_re_5    = S2$RE_5/10000
S2_re_6    = S2$RE_6/10000
S2_re_7    = S2$RE_7/10000
S2_NIR_1   = S2$NIR/10000
S2_NIR_2   = S2$NIR_2/10000
S2_SWIR_1  = S2$SWIR_1/10000
S2_SWIR_2  = S2$SWIR_2/10000


inraster   = stack(NDVI_B8B5, NDVI_B8B6, NDVI_B8B7,NDVI_B8aB5, NDVI_B8aB6, NDVI_B8aB7, PSRI, Clre, NDre1, NDre2,
                 NDWI_1,NDWI_2,NDWI_3, S2_blue, S2_green, S2_red,S2_re_5, S2_re_6, S2_re_7, S2_NIR_1, S2_NIR_2,S2_SWIR_1,S2_SWIR_2 )


names(inraster) = c('NDVI_B8B5', 'NDVI_B8B6', 'NDVI_B8B7','NDVI_B8aB5', 'NDVI_B8aB6', 'NDVI_B8aB7', 'PSRI', 'Clre', 'NDre1', 'NDre2',
              'NDWI_1','NDWI_2','NDWI_3','blue', 'green', 'red', 'RE_5','RE_6', 'RE_7', 'NIR', 'NIR_2', 'SWIR_1', 'SWIR_2')

inraster.rm <- dropLayer(inraster, c('NDVI_B8B7', 'PSRI', 'NDVI_B8B5'))

#inraster.rm = subset(inraster, select=-c('NDVI_B8B7', 'PSRI', 'NDVI_B8B5')) #inraster[-c('NDVI_B8B7', 'PSRI', 'NDVI_B8B5')]
#==================================================================================
# Import training and validation data
#==================================================================================
trainingData   =  shapefile("Training_SamplesPnts.shp")
validationData = shapefile("Validation_SamplesPnts.shp")

#==================================================================================
# Extract raster values for the training samples 
#==================================================================================
training_data  = extract(inraster.rm, trainingData)
training_response = as.factor(trainingData$class) 

#==================================================================================
#Select the number of input variables(i.e. predictors, features)
#==================================================================================
selection<-c(1:20) 
training_predictors = training_data[,selection] 

#==================================================================================
# Train the random forest
#==================================================================================

ntree = 1000    #number of trees to produce per iteration
mtry = 6       # number of variables used as input to split the variables
r_forest = randomForest(training_predictors, y=training_response, ntree = ntree, keep.forest=TRUE, mtry=mtry,importance = TRUE, proximity=TRUE) 

#===================================================================================
#Investigate the OOB (Out-Of-the bag) error
#===================================================================================
r_forest

#===================================================================================
# Assessment of variable importance
#===================================================================================
imp =importance(r_forest)  #for ALL classes individually
imp                        #display importance output in console
varImpPlot(r_forest)
varUsed(r_forest)
importance(r_forest)

#===================================================================================
#Evaluate the impact of the ntree on the accuracy
#===================================================================================
tree_nr_evaluation = data.frame(
             Trees=rep(1:nrow(r_forest$err.rate), times=5), 
                Type=rep(c("OOB", "Corn", "Rice", "Soybean", "Wheat"),
                         each=nrow(r_forest$err.rate)),
               Error = c(r_forest$err.rate[, "OOB"],
                         r_forest$err.rate[, "Corn"],
                          r_forest$err.rate[, "Rice"],
                          r_forest$err.rate[, "Soybean"],
                          r_forest$err.rate[, "Wheat"]))

ggplot(data=tree_nr_evaluation, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))


#=======================================================================================
#Evaluate the impact of the mtry on the accuracy
#========================================================================================
 
 mtry <- tuneRF(training_predictors,training_response, ntreeTry=ntree,
                step=1.5,improve=0.01, trace=TRUE, plot=TRUE, mtryStart = 2)
 best.m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
 print(mtry)
 print(best.m)
 
 #======================================================================================
 #Number of tree nodes
 #======================================================================================
 hist(treesize(r_forest), main= "Number of Nodes for the trees",
      col= "green")

 getTree(r_forest, 1, labelVar = TRUE) # inspect the tree characteristics

#==========================================================================================
# Classify the entire image:define raster data to use for classification
#=========================================================================================
predictor_data = subset(inraster.rm, selection)	
#setwd("......")#change this to your output directory if different 

#==========================================================================================
# Classify the entire image
#=========================================================================================

predictions = predict(predictor_data, r_forest, format=".tif", overwrite=TRUE, progress="text", type="response") 
writeRaster(predictions, 'prediction.tif')
 
#==========================================================================================
# Assess the classification accuracy
#=========================================================================================
validationd = extract(inraster.rm, validationData)
validation=extract(predictions, validationData) # extracts the value of the classified raster at the validation point locations

confusionMatrix(as.factor(validation), as.factor(validationData$classID) )
confusionMatrix(as.factor(validation), as.factor(validationData$classID) )$byClass[, 1]
confusionMatrix(as.factor(validation), as.factor(validationData$classID) )$byClass[]
#==========================================================================================
# Save the classification results
#=========================================================================================

Class_Results = writeRaster(predictions, 'S2_class.tif', overwrite=TRUE) 


#Vizualize samples based on proximity values 
r_forest$proximity
distance.matrix<-dist(1-r_forest$proximity)
distance.matrix
mds.studd<-cmdscale(distance.matrix, eig = TRUE, x.ret = TRUE)
mds.studd
mds.var.per<-round(mds.studd$eig/sum(mds.studd$eig)*100,1)
mds.var.per
mds.values<-data.frame(mds.studd$points)
mds.values
mds.values$Class<-training_response
mds.values$Class
mds.values$Sample<-seq(1,nrow(mds.values),by=1)
mds.values
pplot<-ggplot(data=mds.values,aes(x=X1,y=X2, label=Sample)) + geom_point() + geom_text(aes(color=Class))
pplot 

#MDSplot (r_forest, training_response)

training_predictor

