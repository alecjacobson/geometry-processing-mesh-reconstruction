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
  D.reserve(nx * ny * nz * 2);
  D.setZero();
  
  if (dir == 0) {
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
    D.reserve(nx * ny * nz * 2);
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
