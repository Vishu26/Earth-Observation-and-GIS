#==================================================================================
# R script for practical class Basic Image Processing
# Thresholding, filters
# 
# Last modified: 23 April 2019
#
# This code is developed for educational purpose, 
# to be used in Q4 course "Image Analysis" of MSc course in Geoinformatics at ITC.
# The code is developed by Valentyn Tolpekin, v.a.tolpekin@utwente.nl.
# The code is optimized for clarity rather than for computational efficiency.
#
# Do not remove this announcement.
# The code is distributed "as is", WITH NO WARRANTY whatsoever!
#==================================================================================

#==================================================================================
# Block 1: Open an image
#==================================================================================

rm(list=ls(all=TRUE))
graphics.off()

require(rgdal)
require(Rcpp)

# Location of project data
Root <- getwd()

Path_in <- paste(Root, "/Input", sep="")

setwd(Root)
source("BIP_lib.R")
sourceCpp("BIP_lib.cpp")
library(jpeg)
filename <- "../Google_Earth.jpg"

A <- readJPEG(filename) 

# CHeck image dimensions and subset if too large.
d <- A@grid@cells.dim
nrows0 <- d[1]
ncols0 <- d[2]

# Subset size
nrows <- 500
ncols <- 450

if(nrows < nrows0 & ncols < ncols0) {
	start_row <- round(runif(1, min = 1, max=nrows0 - nrows))
	start_col <- round(runif(1, min = 1, max=ncols0 - ncols))

	end_row <- start_row + nrows
	end_col <- start_col + ncols

	A <- A[start_col:end_col, start_row:end_row]
}

Nb <- dim(A@data)[2]

# Convert the multiband image to grayscale
if(ncol(A@data)>1) A$intensity <- rowMeans(A@data[,]) else names(A@data) <- "intensity"

A.image <- A
# apply histogram strecth
A.image$intensity   <- histstretch(A.image$intensity)

windows()
image(A.image, attr="intensity", col=gray((0:255)/255), axes=TRUE)
title("Intensity image")


# Image dimensions
d <- A@grid@cells.dim
psize <- A@grid@cellsize

nrows <- d[1]
ncols <- d[2]

#==================================================================================
# Block 2: Summary, histogram
#==================================================================================
str(A)
summary(A)

windows()
hist(A$intensity,breaks=100, main = "Histogram")

summary(A$intensity)
#==================================================================================
# Block 3: Thresholding
#==================================================================================
# Manual setting the threshold value
th <- 25
#==================================================================================
# Intermeans
#==================================================================================
# y <- A$intensity

# th <- mean(y)

# eps <- 1.0e-3

# while(TRUE) {
	# condition <- (y<th)

	# y1 <- y[condition]
	# y2 <- y[!condition]

	# m1 <- mean(y1)
	# m2 <- mean(y2)

	# th1 <- mean(c(m1, m2))
	# if(abs(th1-th)<eps) break

	# th <- th1

# }

T <- A

ind <- which(A$intensity >= th)
T$intensity[ind] <- 0
T$intensity[-ind] <- 1

windows()
par(mfrow=c(1,2))
image(T, attr="intensity", col=gray((0:255)/255), axes=TRUE)
title(paste("Thresholded image, th=", th, sep=""))
hist(A$intensity,breaks=100, main = "Histogram")
abline(v=th, lty="dashed", col="blue", lwd=2)

# # Optional
# # Convert thresholded image into a matrix 
# # and save into a file (used later in mathematicl morphology)
# T_mat <- matrix(T@data$intensity,byrow=FALSE,nrow=nrows,ncol=ncols)
# save(T_mat, file=paste(Path_in,"/","thresholded_image.RData", sep=""))

#==================================================================================
# Block 4: Apply linear convolution filtering
#==================================================================================
# Define filter kernel
# Size of the kernel; Only odd numbers
k_nrows <- 5
k_ncols <- 5

# Half width of the kernel, excepth the central pixel
hw_row <- (k_nrows-1)/2
hw_col <- (k_ncols-1)/2

# Matrix of kernel weights
W <- array(0,c(k_nrows,k_ncols))
# Equal kernel values
W[,] <- 1

# Apply gain correction
W <- W /(sum(abs(W)))

# Print the matrix of weights (kernel)
W

# convert image into a matrix
P <- array(0, c(nrows, ncols))
P[,] <- matrix(A@data$intensity,byrow=FALSE,nrow=nrows,ncol=ncols)

# Apply convolutional filter
A$F <- as.vector(convolve_2d_R(P, W))

windows()
par(mfrow=c(2,2))
image(A, attr="intensity", col=gray((0:255)/255), axes=TRUE)
title("Original")
hist(A$intensity, breaks=100, main="Histogram of original")
image(A, attr="F", col=gray((0:255)/255), axes=TRUE)
title("Filtered")
hist(A$F, breaks=100, main="Histogram of filtered")

summary(A$intensity)
summary(A$F)
#==================================================================================
# Block 5: Filtering using Rcpp
#==================================================================================
# A$f <- convolve_2d(A@data$intensity, nrows, ncols, W)
A$f <- as.vector(convolve_2d(P, W))

windows()
par(mfrow=c(2,2))
image(A, attr="intensity", col=gray((0:255)/255), axes=TRUE)
title("Original")
hist(A$intensity, breaks=100, main="Histogram of original")
image(A, attr="f", col=gray((0:255)/255), axes=TRUE)
title("Filtered")
hist(A$f, breaks=100, main="Histogram of filtered")

summary(A$intensity)
summary(A$f)

# Compare results of R and C++ filtering functions
summary(A$F-A$f)
#==================================================================================
# Block 6: Use other filters
#==================================================================================

#==================================================================================
# Manually edit the kernel
#==================================================================================
W[,] <- 1
W <- edit(W)

k_nrows <- dim(W)[1]
k_ncols <- dim(W)[2]

# Half width of the kernel, excepth the central pixel
hw_row <- (k_nrows-1)/2
hw_col <- (k_ncols-1)/2

# Display the kernel as an image
windows()
image(1:k_nrows, 1:k_ncols, W, col=gray((0:255)/255), main="Kernel displayed as an image", xlab="", ylab="")

A$f <- as.vector(convolve_2d(P, W))

windows()
par(mfrow=c(2,2))
image(A, attr="intensity", col=gray((0:255)/255), axes=TRUE)
title("Original")
hist(A$intensity, breaks=100, main="Histogram of original")
image(A, attr="f", col=gray((0:255)/255), axes=TRUE)
title("My filter")
hist(A$f, breaks=100, main="Histogram of filtered")

# # Residual image: filtered - original
# A$residual <- A$f - A$intensity
# windows()
# image(A, attr="residual", col=gray((0:255)/255), axes=TRUE)
# title("Residual image")


#==================================================================================
# Gaussian smoothing
#==================================================================================
# Define filter kernel
# Size of the kernel; Only odd numbers
k_nrows <- 15
k_ncols <- 15

# Half width of the kernel, excepth the central pixel
hw_row <- (k_nrows-1)/2
hw_col <- (k_ncols-1)/2

# Matrix of kernel weights
W <- array(0,c(k_nrows,k_ncols))

# Gaussian
sigma <- (hw_row-1)/6
W[,] <- (row(W)-hw_row-1)^2 + (col(W)-hw_col-1)^2
W <- exp(-W / (2*sigma^2))
W <- W / (2*pi*sigma^2)

W <- W/sum(W)

# Display the kernel as an image
windows()
image(1:k_nrows, 1:k_ncols, W, col=gray((0:255)/255), main="Kernel displayed as an image", xlab="", ylab="")


A$f <- as.vector(convolve_2d(P, W))
windows()
par(mfrow=c(2,2))
image(A, attr="intensity", col=gray((0:255)/255), axes=TRUE)
title("Original")
hist(A$intensity, breaks=100, main="Histogram of original")
image(A, attr="f", col=gray((0:255)/255), axes=TRUE)
title("Filtered")
hist(A$f, breaks=100, main="Histogram of filtered")

# # Residual image: filtered - original
# A$residual <- A$f - A$intensity
# windows()
# image(A, attr="residual", col=gray((0:255)/255), axes=TRUE)
# title("Residual image")

#==================================================================================
# Gaussian gradient
#==================================================================================
# Define filter kernel
# Size of the kernel; Only odd numbers
k_nrows <- 15
k_ncols <- 15

# Half width of the kernel, excepth the central pixel
hw_row <- (k_nrows-1)/2
hw_col <- (k_ncols-1)/2

# Matrix of kernel weights
W <- array(0,c(k_nrows,k_ncols))

# Gaussian gradient
sigma <- (hw_row-1)/6
W[,] <- (row(W)-hw_row-1)^2 + (col(W)-hw_col-1)^2
W <- exp(-W / (2*sigma^2))
W <- W / (2*pi*sigma^2)

W <- -W * (row(W)-hw_row-1) / (sigma^2)

# Display the kernel as an image
windows()
image(1:k_nrows, 1:k_ncols, W, col=gray((0:255)/255), main="Kernel displayed as an image", xlab="", ylab="")


A$f <- as.vector(convolve_2d(P, W))
windows()
par(mfrow=c(2,2))
image(A, attr="intensity", col=gray((0:255)/255), axes=TRUE)
title("Original")
hist(A$intensity, breaks=100, main="Histogram of original")
image(A, attr="f", col=gray((0:255)/255), axes=TRUE)
title("Gaussian gradient x-direction filter")
hist(A$f, breaks=100, main="Histogram of filtered")
#==================================================================================
# Block 7: Non linear filters
#==================================================================================
#==================================================================================
# Robert's filter
#==================================================================================
A$f <- as.vector(roberts(P))
windows()
par(mfrow=c(2,2))
image(A, attr="intensity", col=gray((0:255)/255), axes=TRUE)
title("Original")
hist(A$intensity, breaks=100, main="Histogram of original")
image(A, attr="f", col=gray((0:255)/255), axes=TRUE)
title("Robert's filter")
hist(A$f, breaks=100, main="Histogram of filtered")

# # Residual image: filtered - original
# A$residual <- A$f - A$intensity
# windows()
# image(A, attr="residual", col=gray((0:255)/255), axes=TRUE)
# title("Residual image")

#==================================================================================
# Sobel's filter
#==================================================================================
# Define filter kernel
# Size of the kernel; Only odd numbers
k_nrows <- 3
k_ncols <- 3

# Half width of the kernel, excepth the central pixel
hw_row <- (k_nrows-1)/2
hw_col <- (k_ncols-1)/2

# Matrix of kernel weights
W <- array(0,c(k_nrows,k_ncols))

# Approximate directional first derivatives (Glasbey and Horgan, page 39)
W[1,] <- c(-1, 0, 1)
W[2,] <- c(-2, 0, 2)
W[3,] <- W[1,]

A$grad_x <- as.vector(convolve_2d(P, W))

W[1,] <- c(-1, -2, -1)
W[2,] <- c(0, 0, 0)
W[3,] <- -W[1,]

A$grad_y <- as.vector(convolve_2d(P, W))

A$grad <- sqrt(A$grad_x^2+A$grad_y^2)
A$phi <- atan2(A$grad_y, A$grad_x)

windows()
par(mfrow=c(2,2))
image(A, attr="grad_x", col=gray((0:255)/255), axes=TRUE)
title("Gradient x")
image(A, attr="grad_y", col=gray((0:255)/255), axes=TRUE)
title("Gradient y")
image(A, attr="grad", col=gray((0:255)/255), axes=TRUE)
title("Gradient magnitude, Sobel filter")
image(A, attr="phi", col=gray((0:255)/255), axes=TRUE)
title("Gradient orientation")

#==================================================================================
# Canny filter
#==================================================================================
phi_mat <- array(0, c(nrows, ncols))
phi_mat[,] <- matrix(A$phi, byrow=FALSE, nrow=nrows, ncol=ncols)

A$f <- as.vector(canny(P, phi_mat))
windows()
par(mfrow=c(2,2))
image(A, attr="intensity", col=gray((0:255)/255), axes=TRUE)
title("Original")
hist(A$intensity, breaks=100, main="Histogram of original")
image(A, attr="f", col=gray((0:255)/255), axes=TRUE)
title("Canny filter")
hist(A$f, breaks=100, main="Histogram of filtered")

#==================================================================================
# The END
#==================================================================================
