#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  int ind, i, j, k;
  
  if (dir == 0) {
    D.resize( (nx - 1)*ny*nz, nx*ny*nz );
    for (i = 1; i < nx; i++) {
      for (j = 0; j < ny; j++) {
        for (k = 0; k < nz; k++) {
          ind = i + nx*(j + k * ny);
          D.insert(ind, i-1) = -1 / h;
          D.insert(ind, i) = 1 / h;
        }
      }
    }
  } else if (dir == 1) {
    D.resize( nx*(ny-1)*nz, nx*ny*nz );
    for (i = 0; i < nx; i++) {
      for (j = 1; j < ny; j++) {
        for (k = 0; k < nz; k++) {
          ind = i + nx*(j + k * ny);
          D.insert(ind, j-1) = -1 / h;
          D.insert(ind, j) = 1 / h;
        }
      }
    }
  } else {
    D.resize( nx*ny*(nz-1), nx*ny*nz );
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
        for (k = 1; k < nz; k++) {
          ind = i + nx*(j + k * ny);
          D.insert(ind, k-1) = -1 / h;
          D.insert(ind, k) = 1 / h;
        }
      }
    }
  }
}
