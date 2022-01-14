#==============================================================================
# Binary Morphology Functions for R   	    
# This library was developed for education purposes.
# Algortihms are not optimized for efficiency, but for clarity and simplicity.
# Written by Valentyn Tolpekin		    
# Last modified: December 2017
#==============================================================================

# Global parameters
ShowAxes <- TRUE

# Displays a set in an image
mm_image <- function(A){
	d <- dim(A)
	nrows <- d[1]
	ncols <- d[2]
	
	if(!ShowAxes)par(mai=c(0,0,0,0))
	image(1:nrows,1:ncols,mm_zero(nrows, ncols), col=gray(1:0),xlab="",ylab="",xlim=c(0.5,nrows+0.5),ylim=c(0.5,ncols+0.5))
	rect(0.5,0.5,nrows+0.5,ncols+0.5)

	

	for(i in 1:nrows)
	for(j in 1:ncols)
	{
	
		if((i==1)|(i==nrows)|(j==1)|(j==ncols))
		{
			if(A[i,j])
			{
				lwidth<-1
				lines(c(i-0.5,i-0.5),c(j-0.5,j+0.5),lty=1,lwd=lwidth)
				lines(c(i+0.5,i+0.5),c(j-0.5,j+0.5),lty=1,lwd=lwidth)
				lines(c(i-0.5,i+0.5),c(j-0.5,j-0.5),lty=1,lwd=lwidth)
				lines(c(i-0.5,i+0.5),c(j+0.5,j+0.5),lty=1,lwd=lwidth)		
			}
		}else if(A[i,j]==1)
		{
			if(A[i-1,j]==1)
			{
				lwidth <- 1
			} else lwidth <- 3
			lines(c(i-0.5,i-0.5),c(j-0.5,j+0.5),lty=1,lwd=lwidth)

			if(A[i+1,j]==1)
			{
				lwidth <- 1
			} else lwidth <- 3
			lines(c(i+0.5,i+0.5),c(j-0.5,j+0.5),lty=1,lwd=lwidth)

			if(A[i,j-1]==1)
			{
				lwidth <- 1
			} else lwidth <- 3
			lines(c(i-0.5,i+0.5),c(j-0.5,j-0.5),lty=1,lwd=lwidth)

			if(A[i,j+1]==1)
			{
				lwidth <- 1
			} else lwidth <- 3
			lines(c(i-0.5,i+0.5),c(j+0.5,j+0.5),lty=1,lwd=lwidth)
		}
	}
}


# Displays a SE
mm_display_SE <- function(SE){
	d <- dim(SE)
	nrow_se <- d[1]
	ncol_se <- d[2]
	
	if(!ShowAxes)par(mai=c(0,0,0,0))
	image(1:nrow_se,1:ncol_se,array(0,c(nrow_se,ncol_se)), col=gray(1:0),xlab="",ylab="",xlim=c(0.5,ncol_se+0.5),ylim=c(0.5,ncol_se+0.5))
	rect(0.5,0.5,nrow_se+0.5,ncol_se+0.5)

	for(i in 1:nrow_se)
	for(j in 1:ncol_se) {
		if((i==1)|(i==nrow_se)|(j==1)|(j==ncol_se)) {
			if(!is.na(SE[i,j])){
				if(SE[i,j]==1){
					l_col = "blue"
				} else l_col = "red"
				
				lwidth <- 3
				lines(c(i-0.5,i-0.5), c(j-0.5,j+0.5), lty = 1, lwd = lwidth, col=l_col)
				lines(c(i+0.5,i+0.5), c(j-0.5,j+0.5), lty = 1, lwd = lwidth, col=l_col)
				lines(c(i-0.5,i+0.5), c(j-0.5,j-0.5), lty = 1, lwd = lwidth, col=l_col)
				lines(c(i-0.5,i+0.5), c(j+0.5,j+0.5), lty = 1, lwd = lwidth, col=l_col)
			}
			
		} else {
			if(!is.na(SE[i,j]))
			if(SE[i,j]==1) {
				l_col = "blue"
				if(!is.na(SE[i-1,j])) {
					lwidth <- 1
				} else lwidth <- 3
				lines(c(i-0.5,i-0.5),c(j-0.5,j+0.5),lty=1,lwd=lwidth, col=l_col)

				if(!is.na(SE[i+1,j])) {
					lwidth <- 1
				} else lwidth <- 3
				lines(c(i+0.5,i+0.5),c(j-0.5,j+0.5),lty=1,lwd=lwidth, col=l_col)

				if(!is.na(SE[i,j-1])) {
					lwidth <- 1
				} else lwidth <- 3
				lines(c(i-0.5,i+0.5),c(j-0.5,j-0.5),lty=1,lwd=lwidth, col=l_col)

				if(!is.na(SE[i,j+1])) {
					lwidth <- 1
				} else lwidth <- 3
				lines(c(i-0.5,i+0.5),c(j+0.5,j+0.5),lty=1,lwd=lwidth, col=l_col)
			}

			if(SE[i,j]==0) {
				l_col = "red"
				if(!is.na(SE[i-1,j])) {
					lwidth <- 1
				} else lwidth <- 3
				lines(c(i-0.5,i-0.5),c(j-0.5,j+0.5),lty=1,lwd=lwidth, col=l_col)

				if(!is.na(SE[i+1,j])) {
					lwidth <- 1
				} else lwidth <- 3
				lines(c(i+0.5,i+0.5),c(j-0.5,j+0.5),lty=1,lwd=lwidth, col=l_col)

				if(!is.na(SE[i,j-1])) {
					lwidth <- 1
				} else lwidth <- 3
				lines(c(i-0.5,i+0.5),c(j-0.5,j-0.5),lty=1,lwd=lwidth, col=l_col)

				if(!is.na(SE[i,j+1])) {
					lwidth <- 1
				} else lwidth <- 3
				lines(c(i-0.5,i+0.5),c(j+0.5,j+0.5),lty=1,lwd=lwidth, col=l_col)
			}
		}
	}
}


mm_image_old <- function(A){
	d <- dim(A)
	nrows <- d[1]
	ncols <- d[2]
	
	image(1:nrows,1:ncols,A, col=gray(1:0),xlab="",ylab="",xlim=c(0.5,nrows+0.5),ylim=c(0.5,ncols+0.5))
	
	for(i in 0:(nrows+1)) abline(h=i+0.5,lty=2,lwd=0.5)
	for(i in 0:(ncols+1)) abline(v=i+0.5,lty=2,lwd=0.5)
}

#==============================================================================
# Create sets
#==============================================================================
# Zero object
mm_zero <- function (nrows, ncols) {
	C <- array(0, c(nrows, ncols))
	C
}

# Object of all "1"
mm_one <- function () {
	C <- array(1,c(nrows, ncols))
	return(C)
}

# Object with single non-zero element (x,y)
mm_single <- function (x, y) {
	C <- mm_zero(nrows,ncols)
	C[x,y] <- 1
	return(C)
}

# A square at a position x,y and size 2a+1
mm_square <- function(x, y, a){
	C <- mm_zero(nrows,ncols)
	C[(x - a):((x + a)), ((y - a):(y + a))] <- 1
	C
}

# A disk at a position x,y and radius R
mm_disk <- function (x, y, R) {
	C <- mm_square (x, y, R)
	for (i in -R:R) {
		for (j in -R:R) {
			r <- sqrt(i*i + j*j)
			if (r >= R) C[x+i, y+j] <- 0
		}
	}
	C
}


#==============================================================================
# Define operators
#==============================================================================
# Complementation of object
mm_complementation <- function (A) {
	C <- as.integer(!A)
	dim(C) <- c(nrows, ncols)
	C
}

# Intersection (AND) of sets A and B
mm_intersection <- function (A, B) {
	C <- as.integer(A & B)
	dim(C) <- c(nrows, ncols)
	C
}

# Union (OR) of sets A and B
mm_union <- function (A, B) {
	C <- as.integer(A | B)
	dim(C) <- c(nrows, ncols)
	C
}

# Translation of set A by x,y. 
# Warning: watch the borders!
mm_translation <- function (A, x, y) {
	d <- dim(A)
	nrows <- d[1]
	ncols <- d[2]

	C <- mm_zero(nrows,ncols)
	if ((x>=-nrows+1) & (x<=nrows) & (y>=-ncols+1) & (y<=ncols)) {
		for (i in 1:nrows)
		  for (j in 1:ncols) {
			if ((i+x>=0) & (i+x<=nrows) & (j+y>=0) & (j+y<=ncols)) C[i+x, j+y] <- A[i,j]
		  }
	}
	C
}

# Transposition of structuring element A with respect to point x,y. 
# Warning: no periodicity, watch the borders!
mm_transposition <- function (A, x, y) {
	C <- mm_zero(nrows,ncols)
	for(i in 1:nrows)
	for(j in 1:ncols){
		if((2*x-i>0)&(2*x-i<=nrows)&(2*y-j>0)&(2*y-j<=ncols))
		if(A[2*x-i,2*y-j]==1)
			C[i,j] <- 1
	}
	C
}



# Checks whether B is a subset of A (no displacement allowed!)
mm_issubset <- function(A,B){
	C <- 1

	for(i in 1:nrows)
	for(j in 1:ncols){
		if(B[i,j]==1)
		if(A[i,j]!=1)
			C <- 0
	}
	C
}

# Checks whether B is identical to A (no displacement allowed!)
mm_isidentical <- function (A, B) {

	if(sum((A-B)^2)>0) return(0)
	else return(1)
}

# Erosion of A by B. The reference pixel of B is given by (x,y)
mm_erosion <- function (A, B, x, y) {
	C <- mm_zero(nrows,ncols)
	D <- mm_transposition(B, x, y)

	imin <- 1
	imax <- nrows
	jmin <- 1
	jmax <- ncols

	for (i in 1:nrows) {
		if (sum(D[i,])>0) {
			imin <- i
			break
		}
	}

	for (i in nrows:1) {
		if (sum(D[i,]) > 0) {
			imax <- i
			break
		}
	}

	for (j in 1:ncols) {
		if (sum(D[,j]) > 0) {
			jmin <- j
			break
		}
	}

	for (j in ncols:1) {
		if (sum(D[,j]) > 0) {
			jmax <- j
			break
		}
	}

	if((imax < imin) || (jmax < jmin)) return(C)

	C <- A

	for (i in imin:imax)
	for (j in jmin:jmax) {
		if (D[i,j] > 0) {
			F <- mm_translation(A, i-imin, j-jmin)
			C <- mm_intersection(C, F)
		}
	}

	C <- mm_translation(C, -x+imin, -y+jmin)
	C
}

#  Erosion of the set A by a square with side a
#  the reference point of the square is the top left corner (1,1)
mm_erosionsqr <- function(A,a){

	C <- A

	if(a<=1) return(C)

	for(i in 1:a)
	for(j in 1:a)
	{
		B <- mm_translation(A,1-i,1-j)
		C <- mm_intersection(C,B)

	}

	C <- mm_translation(C,a/2,a/2)

	return(C)
}
##############################################


#
#  Erosion of the set A by a disk with radius R
#  the reference point of the disk is it's center: (R,R)
#

mm_erosiondisk <- function(A,R){

	C <- A

	if(R<=1) return(C)

	D <- mm_disk(R,R,R)

	for(i in 1:(2*R-1))
	for(j in 1:(2*R-1))
	{
		if(D[i,j]>0)
		{
			B <- mm_translation(A,R-i,R-j)
			C <- mm_intersection(C,B)
		}

	}

	return(C)
}
##############################################



#
# Dilation of A by B. The reference pixel of B is given by (x,y)
#


mm_dilation <- function(A,B,x,y){

	C <- mm_zero(nrows,ncols)

	D <- mm_transposition(B,x,y)

	imin<-1
	imax<-nrows
	jmin<-1
	jmax<-ncols

	for(i in 1:nrows)
	{
		if(sum(D[i,])>0)
		{
			imin<-i
			break
		}

	}


	for(i in nrows:1)
	{
		if(sum(D[i,])>0)
		{
			imax<-i
			break
		}

	}

	for(j in 1:ncols)
	{
		if(sum(D[,j])>0)
		{
			jmin<-j
			break
		}

	}


	for(j in ncols:1)
	{
		if(sum(D[,j])>0)
		{
			jmax<-j
			break
		}

	}



	if((imax<imin)||(jmax<jmin)) return(C)

	C <- A

	for(i in imin:imax)
	for(j in jmin:jmax)
	{
		if(D[i,j]>0)
		{
			F <- mm_translation(A, i-imin, j-jmin)
			C <- mm_union(C,F)
		}

	}

	C <- mm_translation(C,-x+imin,-y+jmin)
	return(C)
}
##############################################



#
# Dilation of A by B. The reference pixel of B is given by (x,y)
# Works slower compared to mm_dilation(). Produces artefatcs at the boundaries! To be used only in demo purposes!

mm_dilation2 <- function(A,B,x,y){

	C <- mm_complementation(A)

	C <- mm_erosion(C,B,x,y)

	C <- mm_complementation(C)

	return(C)
}
##############################################



#
# Dilation of A with a square with size a. The reference is (1,1)
#

mm_dilationsqr <- function(A,a){

	C <- A

	if(a<=1) return(C)

	for(i in 1:a)
	for(j in 1:a)
	{
		B <- mm_translation(A,1-i,1-j)
		C <- mm_union(C,B)

	}

	C <- mm_translation(C,a/2,a/2)

	return(C)
}
##############################################


#
#  Dilation of the set A by a disk with radius R
#  the reference point of the disk is it's center: (R,R)
#

mm_dilationdisk <- function(A,R){

	C <- A

	if(R<=1) return(C)

	D <- mm_disk(R,R,R)

	for(i in 1:(2*R-1))
	for(j in 1:(2*R-1))
	{
		if(D[i,j]==1)
		{
			B <- mm_translation(A,R-i,R-j)
			C <- mm_union(C,B)
		}

	}

	return(C)
}
##############################################



#
# Opening of A by B. The reference pixel of B is given by (x,y)
#

mm_opening <- function(A,B,x,y){

	C <- mm_erosion(A,B,x,y)

	D <- mm_transposition(B,x,y)

	C <- mm_dilation(C,D,x,y)

	return(C)
}
##############################################



#
# Opening of A with a square a.
#

mm_openingsqr <- function(A,a){

	C <- mm_erosionsqr(A,a)

	C <- mm_dilationsqr(C,a)

	return(C)
}
##############################################



#
#  Opening of the set A by a disk with radius R
#  the reference point of the disk is it's center: (R,R)
#

mm_openingdisk <- function(A,R){

	C <- mm_erosiondisk(A,R)

	C <- mm_dilationdisk(C,R)

	return(C)
}
##############################################



#
# Closing of A by B. The reference pixel of B is given by (x,y)
#

mm_closing <- function(A,B,x,y){

	C <- mm_dilation(A,B,x,y)

	D <- mm_transposition(B,x,y)

	C <- mm_erosion(C,D,x,y)

	return(C)
}
##############################################


#
# Opening of A with a square of size a.
#

mm_closingsqr <- function(A,a){

	C <- mm_dilationsqr(A,a)

	C <- mm_erosionsqr(C,a)

	return(C)
}
##############################################

#
#  Closing of the set A by a disk with radius R
#  the reference point of the disk is it's center: (R,R)
#

mm_closingdisk <- function(A,R){

	C <- mm_dilationdisk(A,R)

	C <- mm_erosiondisk(C,R)

	return(C)
}
##############################################



#
#  Erosion of the set A by a set of two points separated by 
#  distance a with relative orientation given by parameters tx,ty 
#  tx = 1, ty = 0:   Horizontally
#  tx = 0, ty = 1:   Vertically
#  tx = 1, ty = 1:   Diagonal (NE-SW)
#  tx = 1, ty = -1:  Diagonal (SE-NW)
#

mm_erosiontxtr <- function(A,a,tx,ty){

	C <- mm_zero(nrows, ncols)

	if(a<1) return(C)

	for(i in (a+1):(nrows-a-1))
	for(j in (a+1):(ncols-a-1))
	{
		if(A[i,j]==1)
		{
			if(A[i+tx*a,j+ty*a]==1)
				C[i,j] <- 1
		}

	}

	return(C)
}
#==============================================================================
# Hit and miss transform of A with a structuring element SE
# SE is assumed with odd dimensions
#==============================================================================
mm_hit_and_miss <- function(A, SE) {
	d <- dim(A)
	nrows <- d[1]
	ncols <- d[2]

	d <- dim(SE)
	nrow_se <- d[1]
	ncol_se <- d[2]
	hw_row <- (nrow_se - 1)/2
	hw_col <- (ncol_se - 1)/2

	sev <- as.vector(SE)
	ind <- which(sev!=0)
	sev <- sev[ind]

	C <- mm_zero(nrows, ncols)

	for (i_row in (hw_row+1):(nrows-hw_row)) {
		for (j_col in (hw_col+1):(ncols-hw_col)) {
			sub_a <- as.vector(A[(i_row-hw_row):(i_row+hw_row), (j_col-hw_col):(j_col+hw_col)])
			sub_a <- sub_a[ind]
			sub_a[sub_a==0] <- -1
			C[i_row, j_col] <- identical(sub_a, sev)
		}
	}
	C
}

# Rotation of SE by angle that is integer times 90 degrees
mm_rotate_SE <- function (SE, angle) {
	d <- dim(SE)
	nrow_se <- d[1]
	ncol_se <- d[2]
	hw_row <- (nrow_se - 1)/2
	hw_col <- (ncol_se - 1)/2

	C <- array(0, d)
	C[, ] <- 0

	alpha <- angle/180*pi

	x <- col(SE) - hw_row - 1
	y <- row(SE) - hw_col - 1

	r <- sqrt(x^2 + y^2)
	phi <-  atan2(y, x)

	phi <- phi - alpha
	x1 <- round(r * cos(phi))
	y1 <- round(r * sin(phi))

	x1 <- x1 + 1 + hw_row
	y1 <- y1 + 1 + hw_col
	pn1 <- (as.vector(y1) - 1) * ncol_se + as.vector(x1)
	
	pn <- (as.vector(y+hw_col+1) - 1) * ncol_se + as.vector(x+hw_row+1)
	C[pn1] <- SE[pn]

	C
}

mm_thinning <- function(A) {

	windows()
	mm_image(A)

	windows()
	par(mfrow = c(1, 2))
	C <- A
	SE <- rbind(c(0, 0, 0),
		c(NA,  1, NA),
		c(1, 1, 1))
	B <- mm_hit_and_miss(C, SE)
	C <- C - B
	mm_display_SE(SE)
	mm_image(C)

	SE <- mm_rotate_SE(SE, 90)
	B <- mm_hit_and_miss(C, SE)
	C <- C - B
	mm_display_SE(SE)
	mm_image(C)

	SE <- mm_rotate_SE(SE, 90)
	B <- mm_hit_and_miss(C, SE)
	C <- C - B
	mm_display_SE(SE)
	mm_image(C)

	SE <- mm_rotate_SE(SE, 90)
	B <- mm_hit_and_miss(C, SE)
	C <- C - B
	mm_display_SE(SE)
	mm_image(C)

	windows()
	mm_image(C)


}




