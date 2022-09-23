#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  if (dir == 0) {
	  for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
         	   for (int i = 0; i < nx - 1; i++) {
	 		D.insert(i + (nx - 1) * (j + k * ny), (i+1) + nx * (j + k * ny)) = 1;
	 		D.insert(i + (nx - 1) * (j + k * ny), i + nx * (j + k * ny)) = -1;
        	} 
       	   }
         }
         D = D / h; 
   }

   if (dir == 1) {
         for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny - 1; j++) {
				D.insert(i + nx * (j + k * (ny -1)), i + nx * ((j+1) + k * ny)) = 1; 
				D.insert(i + nx * (j + k * (ny -1)), i + nx * (j + k * ny)) = -1; 
		}
             }
         }
         D = D / h; 
   }

   if (dir == 2) {
	for (int k = 0; k < nz - 1; k++) {
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				D.insert(i + nx * (j + k * ny), i + nx * (j + (k+1) *ny)) = 1;
				D.insert(i + nx * (j + k * ny), i + nx * (j + k * ny)) = -1;
			}
		}
        }
         D = D / h; 
   }
}
