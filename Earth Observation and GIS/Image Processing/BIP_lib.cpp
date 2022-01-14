#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix convolve_2d(NumericMatrix A, NumericMatrix K) {
   
	int n_row = A.nrow(), k_row = K.nrow();
	int n_col = A.ncol(), k_col = K.ncol();
	
	int hw_row = (k_row-1)/2, hw_col = (k_col-1)/2;
	int nr, nc;

	NumericMatrix output(n_row, n_col);

	for(int row = 0; row < n_row; row++) {
		for(int col = 0; col < n_col; col++) {
			for(int drow = -hw_row; drow <= hw_row; drow++) {
				for(int dcol = -hw_col; dcol <= hw_col; dcol++) {
					nr = row-drow;
					nc = col-dcol;
					if(nr>=0 && nr<n_row && nc>=0 && nc <n_col) {
						output(row, col) += A(nr, nc) * K(drow+hw_row, dcol+hw_col);
					}
				}
			}
		}
   }

   return output;
}

// [[Rcpp::export]]
NumericMatrix roberts(NumericMatrix A) {
   
	int n_row = A.nrow();
	int n_col = A.ncol();
	
	NumericMatrix output(n_row, n_col);

	for(int row = 0; row < n_row-1; row++){
		for(int col = 0; col < n_col-1; col++){
			output(row, col) = abs(A(row, col) - A(row+1, col+1)) + abs(A(row+1,col) - A(row, col+1));
		}
   }
   
   return output;
}

// [[Rcpp::export]]
NumericMatrix canny(NumericMatrix A, NumericMatrix phi_mat) {
   
	int n_row = A.nrow();
	int n_col = A.ncol();
	int nr, nc;

	NumericMatrix W(3,3);
	NumericMatrix output(n_row, n_col);
	double phi;

	for(int row = 0; row < n_row-1; row++){
		for(int col = 0; col < n_col-1; col++){
			phi = phi_mat(row, col);
			W(0, 0) = 0.5*sin(phi);
			W(0, 1) = pow(sin(phi), 2.0);
			W(0, 2) = -0.5*sin(phi)*cos(phi);
			W(1, 0) = pow(cos(phi), 2.0);
			W(1, 1) = -2.0;
			W(1, 2) = pow(cos(phi), 2.0);
			W(2, 0) = -0.5*sin(phi)*cos(phi);
			W(2, 1) = pow(sin(phi), 2.0);
			W(2, 2) = 0.5*sin(phi)*cos(phi);

			for(int drow = -1; drow <= 1; drow++) {
				for(int dcol = -1; dcol <= 1; dcol++) {
					nr = row-drow;
					nc = col-dcol;
					if(nr>=0 && nr<n_row && nc>=0 && nc <n_col) {
						output(row, col) = A(nr, nc) * W(drow+1, dcol+1);
					}
				}
			}
			
		}
   }
   
   return output;
}

