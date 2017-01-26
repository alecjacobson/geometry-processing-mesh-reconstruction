#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  int row_ind, col_ind1, col_ind2, i, j, k;
  D.reserve(nx * ny * nz * 2);
  D.setZero();
  
  if (dir == 0) {
    for (i = 0; i < nx-1; i++) {
      for (j = 0; j < ny; j++) {
        for (k = 0; k < nz; k++) {
          row_ind = i + (nx - 1)*(j + k * ny);
          col_ind1 = i + nx * (j + k * ny);
          col_ind2 = (i + 1) + nx * (j + k * ny);
          D.insert(row_ind, col_ind1) = -1 / h;
          D.insert(row_ind, col_ind2) = 1 / h;
        }
      }
    }
  } else if (dir == 1) {
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny-1; j++) {
        for (k = 0; k < nz; k++) {
          row_ind = i + nx*(j + k * (ny-1));
          col_ind1 = i + nx*(j + k * ny);
          col_ind2 = i + nx*( (j+1) + k * ny);
          D.insert(row_ind, col_ind1) = -1 / h;
          D.insert(row_ind, col_ind2) = 1 / h;
        }
      }
    }
  } else {
    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
        for (k = 0; k < nz-1; k++) {
          row_ind = i + nx*(j + k * ny);
          col_ind1 = i + nx*(j + k * ny);
          col_ind2 = i + nx*(j + (k+1) * ny);
          D.insert(row_ind, col_ind1) = -1 / h;
          D.insert(row_ind, col_ind2) = 1 / h;
        }
      }
    }
  }
}
