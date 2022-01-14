#==================================================================================
# R library for practical class Basic Image Processing
# Last modified: 25 April 2019
#
# This code is developed for educational purpose, 
# to be used in Q4 course "Image Analysis" of MSc course in Geoinformatics at ITC.
# The code is developed by Valentyn Tolpekin, v.a.tolpekin@utwente.nl.
# The code is optimized for clarity rather than for computational efficiency.
#
# Do not remove this announcement.
# The code is distributed "as is", WITH NO WARRANTY whatsoever!
#==================================================================================
histstretch <- function(x) {
  cur.lim <- quantile(x, c(0.025,0.975), na.rm=TRUE)
  val <- pmax(cur.lim[1], pmin(cur.lim[2], x))
  val <- floor(255*(val - cur.lim[1])/(cur.lim[2] - cur.lim[1]))
  val[is.na(val)] <- 0
  return(val)
}

convolve_2d_R <- function(P, W) {
	F <- array(0, dim(P))
	nrows <- dim(P)[1]
	ncols <- dim(P)[2]
	k_nrows <- dim(W)[1]
	k_ncols <- dim(W)[2]
	i0 <- (k_nrows-1)/2
	j0 <- (k_ncols-1)/2

	# Dilate the image with zero rows and columns
	Pdil <- array(0, c(nrows+k_nrows-1, ncols+k_ncols-1))
	Pdil[i0 + 1:nrows, j0 + 1:ncols] <- P

	for(i in 1:nrows)
	for(j in 1:ncols) {
		# observe the flipping of the image wrt kernel
		F[i,j] <- sum(Pdil[i0+(i+i0):(i-i0), j0+(j+j0):(j-j0)] * W)
	}

	# return(as.vector(F))
	return(F)
}
