#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix neighbourhood_systemC(int, int);
NumericVector Grow_region_seed(NumericVector, NumericVector, int, int, int, double, double);

// This function computes hit and miss transform of a binary image A
// by a strctural element se, values 1 (hit),0 (don't care),-1 (miss)
// se dimensions are odd numbers

// [[Rcpp::export]]
IntegerMatrix hit_and_miss(IntegerMatrix A, IntegerMatrix se){

	int n_row = A.nrow(), n_col = A.ncol();
	int se_nrow = se.nrow(), se_ncol = se.ncol();
	int mid_row = (se_nrow - 1) / 2;
	int mid_col = (se_ncol - 1) / 2;
	int a, b, count;

	IntegerMatrix out(n_row, n_col);

	for (int row = mid_row; row < n_row - mid_row; row++) {
		for (int col = mid_col; col < n_col-mid_col; col++) {
			count = 0;

			for (int i_row = -mid_row; i_row <= mid_row; i_row++) {
				for (int i_col = - mid_col; i_col <= mid_col; i_col++) {
					b = se(i_row+mid_row, i_col+mid_col);
					if(b==0) continue;
					a = A(row + i_row, col + i_col);

					if (a==1) {
						if (b==-1) count++;
					}
					if (a==0) {
						if (b==1) count++;
					}
					//Rcout << "row=" << row << "_col=" << col << "_ir=" << i_row << "_ic=" << i_col << "_a=" << a << "_b=" << b << "count=" << count << "\n";
				}
			}

			if(count==0) {
				out(row, col) = 1;
			} else out(row, col) = 0;
		}
	}
   
   return out;
}


// [[Rcpp::export]]
NumericMatrix neighbourhood_systemC(int M, int N){
   
	int npix = M*N;
	int i, j, pn;
	
	NumericMatrix out(npix,9);
   
	for(i = 1; i < M-1; i++){
		for(j = 1; j < N-1; j++){
			pn = j*M + i;
			out(pn,0) = 8;
			out(pn,1) = pn - M - 1;
			out(pn,2) = pn - M;
			out(pn,3) = pn - M + 1;
			out(pn,4) = pn - 1;
			out(pn,5) = pn + 1;
			out(pn,6) = pn + M - 1;
			out(pn,7) = pn + M;
			out(pn,8) = pn + M + 1;
		}
		
		j = 0;
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - 1;
		out(pn,2) = pn + 1;
		out(pn,3) = pn + M - 1;
		out(pn,4) = pn + M;
		out(pn,5) = pn + M + 1;
		
		j = N-1;
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - M - 1;
		out(pn,2) = pn - M;
		out(pn,3) = pn - M + 1;
		out(pn,4) = pn - 1;
		out(pn,5) = pn + 1;
	}
	
	i=0;
	for(j = 1; j < N-1; j++){
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - M;
		out(pn,2) = pn - M + 1;
		out(pn,3) = pn + 1;
		out(pn,4) = pn + M;
		out(pn,5) = pn + M + 1;
	}

	j = 0;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn + 1;
	out(pn,2) = pn + M;
	out(pn,3) = pn + M + 1;
	
	
	j = N - 1;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn - M;
	out(pn,2) = pn - M + 1;
	out(pn,3) = pn + 1;

	i = M - 1;
	for(j = 1; j < N-1; j++){
		pn = j*M + i;
		out(pn,0) = 5;
		out(pn,1) = pn - M - 1;
		out(pn,2) = pn - M;
		out(pn,3) = pn - 1;
		out(pn,4) = pn + M - 1;
		out(pn,5) = pn + M;
	}
	
	j=0;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn - 1;
	out(pn,2) = pn + M - 1;
	out(pn,3) = pn + M;

	j=N-1;
	pn = j*M + i;
	out(pn,0) = 3;
	out(pn,1) = pn - M - 1;
	out(pn,2) = pn - M;
	out(pn,3) = pn - 1;

	return out;
}

// [[Rcpp::export]]
NumericVector Grow_region_seed(NumericVector f_magn, NumericVector f_angle, int M, int N, int seed, double tau, double magn_min){
   
	int npix = M*N;
	int i, j, k, pn, pn2, grow;
	int count=0;
	double alpha, theta, sx, sy;
	NumericVector marked(npix);
	NumericVector seg(npix);
	NumericMatrix narr(npix,9);
   
	narr = neighbourhood_systemC(M,N);
	
	pn = seed-1; // in C++ array counter starts with zero

	seg[count] = pn;
	marked[pn] = 1;

	count++;

	grow = 1;
	
	while(grow==1){
		for(i=0; i<count; i++){
			grow=0;
			pn = seg[i];

			// mean angle of the current segment
			sx = 0.0;
			sy = 0.0;
			for(k=0;k<count;k++){
				sx+=cos(f_angle[seg[k]]);
				sy+=sin(f_angle[seg[k]]);
			}
			theta = atan2(sy,sx);
			
			for(j=1; j<narr(pn,0)+1;j++){
				pn2 = narr(pn,j);
				
				
				alpha = abs(theta - f_angle[pn2]);
			
				alpha += (alpha>PI) ? -2.0*PI : (alpha<-PI) ? (2.0*PI) : 0.0;	
			
				//if(abs(alpha)<=tau && marked[pn2]==0){
				if(abs(alpha)<=tau && f_magn[pn2]>0.0){
					marked[pn2] = 1;
					f_magn[pn2] = 0.0;
					seg[count]=pn2;
					count++;
					grow=1;
				}
			}
		}
	}

	NumericVector out(count+1);

	out[0] = count;
	for(i=0; i<count; i++){
		out[i+1] = seg[i];//+1; // In R array counter starts with 1
	}
	
	return out;
}

// [[Rcpp::export]]
NumericMatrix segm_map(NumericVector f_magn, NumericVector f_angle, int M, int N, double tau, double magn_min, double min_area){
   
	int npix = M*N;
	int i,j,k, nj, pn, pn2, grow;
	int n_left = 0;
	int count = 0;
	double alpha, theta, sx, sy;
	double magn_max = 0.0;
	NumericVector magn(npix);

	NumericMatrix output(npix,2);

	int seg_length, seg_count = 0;
	NumericVector seg(npix), seg_id(npix);
	NumericVector marked(npix);

	NumericMatrix narr(npix,9);
	
	narr = neighbourhood_systemC(M,N);
	
	// mark the low gradient magnitude pixels as inappropriate for region growing
	for(i=0; i<npix; i++){
		magn[i] = f_magn[i];
		if(magn[i]<=magn_min){
			magn[i] = 0;
		}else{
			n_left++;
			if(magn[i]>magn_max) magn_max = magn[i];
		}
	}
	
	if(magn_max<magn_min){
		return 0;
	}
	
	while(n_left>0){
		//seg_count++;
		count = 0;
		
		magn_max = 0.0;
		
		for(i=0; i<npix; i++){
			if(magn[i]>magn_max){
				magn_max = magn[i];
				pn=i;
			}
		}

		seg[count] = pn;
		marked[pn] = 1;
		magn[pn] = 0.0;

		count++;

		grow = 1;

		while(grow==1){
			// mean angle of the current segment
			sx = 0.0;
			sy = 0.0;
			for(k=0;k<count;k++){
				sx+=cos(f_angle[seg[k]]);
				sy+=sin(f_angle[seg[k]]);
			}

			theta = atan2(sy,sx);

			for(i=0; i<count; i++){
				grow=0;
				pn = seg[i];

				for(j=1; j<narr(pn,0)+1;j++){
					pn2 = narr(pn,j);
					if(magn[pn2]>0.0){
						alpha = abs(theta - f_angle[pn2]);
						alpha += (alpha>PI) ? -2.0*PI : (alpha<-PI) ? (2.0*PI) : 0.0;	

						//if(abs(alpha)<=tau && marked[pn2]==0){
						if(abs(alpha)<=tau){
							marked[pn2] = 1;
							magn[pn2] = 0.0;
							seg[count]=pn2;
							count++;
							grow=1;
						}
					}
				}
				
			}
			
		}

		if(count>min_area){
			seg_count++;				
			for(k = 0; k < count; k++){
				nj = seg[k];
				output(nj,0) = seg_count;  // segment id
				output(nj,1) = count;	  // segment size
			}
		}

		// update n_left
		n_left = 0;
		for(i=0; i<npix; i++) if(magn[i]>=magn_min) n_left++;
	}

	return output;
}

