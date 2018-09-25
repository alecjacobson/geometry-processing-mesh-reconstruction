#include "fd_partial_derivative.h"
#include <assert.h>

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

    int nx_D = nx, ny_D = ny, nz_D = nz;

    switch (dir) {
    case 1:
        nx_D--;
        break;
    case 2:
        ny_D--;
        break;
    case 3:
        nz_D--;
        break;
    default:
        assert(0);
    }

    D.resize(nx_D * ny_D * nz_D, nx * ny * nz);

    for (int i = 0; i < nx_D; i++) {
        for (int j = 0; j < ny_D; j++) {
            for (int k = 0; k < nz_D; k++) {

                int row_index = i + j * nx_D + k * ny_D * nx_D;
                int col_index_prev = i + j * nx + k * ny * nx;
                int col_index_next;

                switch (dir) {
                case 1:
                    col_index_next = (i + 1) + j * nx + k * ny * nx;
                    break;
                case 2:
                    col_index_next = i + (j + 1) * nx + k * ny * nx;
                    break;
                case 3:
                    col_index_next = i + j * nx + (k + 1) * ny * nx;
                    break;
                }

                D(row_index, col_index_prev) = -1;
                D(row_index, col_index_next) = 1;
            }
        }
    }
}
