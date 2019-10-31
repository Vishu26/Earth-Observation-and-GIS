library(raster)
data = brick('IndianPines.tif')
gt = raster('IndianPinesGt.tif')
img = as.array(data)
X = matrix(img, nrow = data@nrows*data@ncols, ncol = data@data@nlayers)

image(img[,,3])
#image(gt)

c = cov(X)
e = eigen(c)
lamb = e$values
U = e$vectors

A = t(U)

Y = A %*% t(X)

prin_1 = matrix(Y[1,], nrow = data@nrows, ncol = data@ncols)
image(prin_1)
sigma = cov(t(Y))
