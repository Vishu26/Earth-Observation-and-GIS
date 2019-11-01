library(raster)
library(caret)
data = brick('IndianPines.tif')
gt = raster('IndianPinesGt.tif')
img = as.array(data)
X = matrix(img, nrow = data@nrows*data@ncols, ncol = data@data@nlayers)
imggt = as.matrix(gt)
#image(img[,,3])
#image(imggt)

c = cov(X)
e = eigen(c)
lamb = e$values
U = e$vectors

A = t(U)

Y = A %*% t(X)

prin_1 = matrix(Y[1,], nrow = data@nrows, ncol = data@ncols)
#image(prin_1)
#sigma = cov(t(Y))

x = as.matrix(Y[1:5,])
xx = t(x)

tr = as.factor(imggt)
dataset = scale(data.frame(xx))

control <- trainControl(method="cv", number=3)
metric <- "Accuracy"

fit.cart <- train(dataset, tr, method='lda', metric=metric)
print(fit.cart)
