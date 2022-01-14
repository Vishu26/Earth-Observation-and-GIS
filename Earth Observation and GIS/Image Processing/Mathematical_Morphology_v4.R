#==================================================================================
# R script for practical class Basic Image Processing
# Mathematical Morphology: binary images
#
# Last modified: 26 April 2019
#
# This code is developed for educational purpose, 
# to be used in Q4 course "Image Analysis" of MSc course in Geoinformatics at ITC.
# The code is developed by Valentyn Tolpekin, v.a.tolpekin@utwente.nl.
# The code is optimized for clarity rather than for computational efficiency.
#
# Do not remove this announcement.
# The code is distributed "as is", WITH NO WARRANTY whatsoever!
#==================================================================================

rm(list=ls(all=TRUE))
graphics.off()

require(Rcpp)
# Root <- "C:/Programming/M8"
Root <- getwd()

setwd(Root)
source("Binary_Morphology.R")
sourceCpp("binary_mm.cpp")
#==================================================================================
# Manually define the binary image
#==================================================================================
nrows <- 30
ncols <- 30

A <- array(0,c(nrows, ncols))

A[16:24,16:24] <- 1

A[24,19:21] <- 0

A[20,25:29] <- 1
A[18:22,27] <-1
A[18:19,28] <-1

A[19:21,19:20]<-0
A[14:16,19:21] <- 1

A[16:17,24] <- 0
A[19,20] <- 1

A <- mm_translation(A,-4,-5)

A <- A[-(1:5),-(1:5)]

A <- A[-(22:25),-(22:25)]

#==================================================================================
# Manual editing of set A
#==================================================================================
# A <- edit(A)
#==================================================================================
# Optional: read a binary image form previous exercise
#==================================================================================
# Path_in <- paste(Root, "/", "Input", sep="")
# load(file=paste(Path_in, "/", "thresholded_image.RData", sep=""))
# A <- T_mat
#==================================================================================
# End of Optional part
#==================================================================================

nrows <- dim(A)[1]
ncols <- dim(A)[2]

windows()
#par(mfrow=c(2,3))
mm_image(A)
#title("Set A")

#==================================================================================
# Erosion
#==================================================================================
E <- mm_erosionsqr(A,3)

windows()
mm_image(E)
title("Erosion of A by a 3x3 square")

#==================================================================================
# Dilation
#==================================================================================
D <- mm_dilationsqr(A,3)

windows()
mm_image(D)
title("Dilation of A by a 3x3 square")

#==================================================================================
# Opening
#==================================================================================
O <- mm_openingsqr(A,3)

windows()
mm_image(O)
title("Opening of A by a 3x3 square")

#==================================================================================
# Closing
#==================================================================================
C <- mm_closingsqr(A,3)

windows()
mm_image(C)
title("Closing of A by a 3x3 square")


#==================================================================================
# Hit and Miss
#==================================================================================
SE <- rbind(c(0, 1, 0),
		c(1,  1, -1),
		c(0, -1, -1))

C <- mm_hit_and_miss(A, SE)

SE <- rbind(c(0, 1, 0),
		c(-1,  1, 1),
		c(-1, -1, 0))

C2 <- mm_hit_and_miss(A, SE)

C <- mm_union(C, C2)

SE <- rbind(c(0, -1, -1),
		c(1,  1, -1),
		c(0, 1, 0))

C2 <- mm_hit_and_miss(A, SE)

C <- mm_union(C, C2)

SE <- rbind(c(-1, -1, 0),
		c(-1,  1, 1),
		c(0, 1, 0))

C2 <- mm_hit_and_miss(A, SE)

C <- mm_union(C, C2)

windows()
mm_image(C)
title("Hit and miss of A by a 3x3 SE")

#==================================================================================
# Thinning
#==================================================================================
A0 <- A
iter <- 1
windows(record=TRUE)

while (TRUE) {
	C <- A

	fid <- 1
	
	SE <- rbind(c(-1, -1, -1),
		c(0,  1, 0),
		c(1, 1, 1))

	B <- mm_hit_and_miss(C, SE)
	# mm_image(B)

	SE1 <- SE
	SE1[SE1==0] <- -1
	SE1[is.na(SE1)] <- 0
	B1 <- hit_and_miss(C, SE)

	C <- C - B
	# mm_image(C)
	# title(main=paste("id=", fid, sep=""))

	for (fid in 2:4) {
		SE <- mm_rotate_SE(SE, 90)
		B <- mm_hit_and_miss(C, SE)
		C <- C - B
		# mm_image(C)
		# title(main=paste("id=", fid, sep=""))
	}

	fid <- 1

	SE <- rbind(c(0, -1, -1),
		c(1,  1, -1),
		c(0, 1, 0))

	B <- mm_hit_and_miss(C, SE)
	C <- C - B
	# mm_image(C)
	# title(main=paste("id=", fid, sep=""))

	for (fid in 2:4) {
		SE <- mm_rotate_SE(SE, 90)
		B <- mm_hit_and_miss(C, SE)
		C <- C - B
		# mm_image(C)
		# title(main=paste("id=", fid, sep=""))
	}

	A <- C

	if(identical(A, A0)) break

	mm_image(A)
	title(main=paste("iteration=", iter, sep=""))
	iter <- iter + 1

	A0 <- A
}
#==================================================================================
# The End
#==================================================================================
