library(raster)
library(readxl)
library(dtw)
library(ggplot2)
setwd("~/Downloads/dtw_shastry")
strname <- "anand_subset_aau_tobacco_vector_tiff.tif"
fileR <- brick(strname)
mask <- raster("anand_subset_aau_tobacco_vector_tiff_cls.tif")
#print(fileR)
#plot(mask)

pure <- read_excel("pure_profile_tobacco.xlsx")
no <- as.double(unlist(pure["Control_Pb"]))
ppm_5 <- as.double(unlist(pure["5_ppm_Pb"]))
ppm_10 <- as.double(unlist(pure["10_ppm_Pb"]))
ppm_15 <- as.double(unlist(pure["15_ppm_Pb"]))

cla = matrix(0L, nrow=1201, ncol=699)

#query = as.double(unlist(fileR[450, 545]))
#reference = as.double(unlist(ppm_5))
#a = dtw(query, reference)
#print(a["distance"])

for(i in 400:500){
  if (i%%10==0){
    print(i)
  }
  for(j in 400:500){
    if (mask[i, j]==1){
      query = as.double(unlist(fileR[i, j]))
      ind = which.min(c(dtw(query, no)["distance"], dtw(query, ppm_5)["distance"], dtw(query, ppm_10)["distance"], dtw(query, ppm_15)["distance"]))
      cla[i, j] = ind 
    }
  }
}
