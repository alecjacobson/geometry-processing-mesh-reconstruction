#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <igl/cat.h>

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
    std::vector< Eigen::SparseMatrix<double> > D(4);
    for (int i = 0; i < 3; i++) {
        fd_partial_derivative(nx, ny, nz, h, i, D[i]);
    }
    igl::cat(1, D[0], D[1], D[3]);
    igl::cat(1, D[3], D[2], G);
}
