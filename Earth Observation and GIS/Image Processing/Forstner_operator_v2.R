#==================================================================================
# R script for practical class Basic Image Processing
# Forstner interest operator
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

#=============================================================================
# Block 1: variable definitions, data import, preparation
#=============================================================================
rm(list=ls(all=TRUE))
graphics.off()

require(rgdal)
require(Rcpp)

# Show intemrediate results? (slower)
ShowIntermediate <- FALSE

# Filter kernel: n = size
n <- 7
# Threshold value for q: q >= tq, where tq is the quant_q quantile
quant_q <- 0.9
# Threshold value for w: w >= tw, where tw is the quant_w quantile 
quant_w <- 0.9

# Location of root folder and input data
Root <- getwd()
Path_in <- paste(Root, "/Input",sep="")

setwd(Root)
sourceCpp("BIP_lib.cpp")

histstretch <- function(x) {
  cur.lim <- quantile(x, c(0.025,0.975), na.rm=TRUE)
  val <- pmax(cur.lim[1], pmin(cur.lim[2], x))
  val <- floor(255*(val - cur.lim[1])/(cur.lim[2] - cur.lim[1]))
  val[is.na(val)] <- 0
  return(val)
}

filename <- "ITAQUI009013NeighOrtho7700X_080829.jpg"

A <- readGDAL(paste(Path_in, "/", filename, sep="")) 

# Image is too large. Subset.
d <- A@grid@cells.dim
nrows0 <- d[1]
ncols0 <- d[2]

nrows <- 300
ncols <- 300

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

# Apply histogram strecth
A.image@data$intensity   <- histstretch(A.image@data$intensity)

if(ShowIntermediate) {
	windows()
	image(A.image, attr="intensity", col=gray((0:255)/255), axes=TRUE)
	title("Intensity image")
}

# Image dimensions
d <- A@grid@cells.dim
psize <- A@grid@cellsize

M <- d[1]
N <- d[2]
npix <- prod(d)
#=============================================================================
# Compute gradients
#=============================================================================
G <- array(0, c(M, N))
G[,] <- A$intensity

W <- (n-1)/2

# Add zero rows and columns
# Offsetting the image array by 2*w zero values 

# Compute gradients, Gx and Gy
Gxshift <- G
Gxshift[2:M,] <- Gxshift[1:(M-1),]
Gxshift[1,] <- Gxshift[2,]

Gx <- Gxshift-G

Gyshift <- G
Gyshift[,2:N] <- Gyshift[,1:(N-1)]
Gyshift[,1] <- Gyshift[,2]

Gy <- Gyshift-G
#=============================================================================
# Apply spatial averaging of the gradients, box filter
#=============================================================================
# Define matrices N_ij
N11  <- N12  <- N21  <- N22  <- array(0, npix)

# Apply spatial averaging
# window size is n x n
W <- array(1, c(n,n))
W <- W/sum(W)

N11 <- as.vector(convolve_2d(Gx^2, W))
N12 <- as.vector(convolve_2d(Gx*Gy, W))
N21 <- N12
N22 <- as.vector(convolve_2d(Gy^2, W))

# Display the components of matrix N (per pixel)
A$N11 <- N11
A$N12 <- N12
A$N22 <- N22

if(ShowIntermediate) {
	windows()
	par(mfrow=c(2,2))
	image(A,attr="N11", col=gray((0:255)/255),axes=TRUE)
	title(main="N11")
	image(A,attr="N12", col=gray((0:255)/255),axes=TRUE)
	title(main="N12")
	image(A,attr="N22", col=gray((0:255)/255),axes=TRUE)
	title(main="N22")
}

# Compute w and q
detN <- N11*N22-N12*N21
TN <- N11+N22

q <- 4 * detN / TN^2
w <- detN / TN

A$q <- q
A$w <- w

tq <- quantile(q, quant_q, na.rm = TRUE)
tw <- quantile(w, quant_w, na.rm = TRUE)

if(ShowIntermediate) {
	windows()
	par(mfrow=c(2,2))
	image(A, attr="q", col=gray((0:255)/255), axes=TRUE)
	title("q")
	hist(q)
	abline(v=tq, col="blue", lty="dotted")
	image(A, attr="w", col=gray((0:255)/255), axes=TRUE)
	title("w")
	hist(w)
	abline(v=tw, col="blue", lty="dotted")
}

# Apply thresholding on q and w
# Combination of q>tq and w>tw yields distinct point of interest operator
ind <- which(q>tq & w>tw)

windows()
image(A.image, attr="intensity", col=gray((0:255)/255), axes=TRUE)
title(main="Forstner operator", sub=paste("n=",n,", q>=",format(tq,digits=3),", w>=",format(tw,digits=3)," (upper ",100*(1-quant_w),"%)",sep=""))
points(coordinates(A)[ind,],col="red",pch=16,cex=0.5)
#=============================================================================
# The END
#=============================================================================
