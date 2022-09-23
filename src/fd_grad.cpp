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
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

    Eigen::SparseMatrix<double> Dx;
    Eigen::SparseMatrix<double> Dy;
    Eigen::SparseMatrix<double> Dz;

    fd_partial_derivative(nx, ny, nz, h, 1, Dx);
    fd_partial_derivative(nx, ny, nz, h, 2, Dy);
    fd_partial_derivative(nx, ny, nz, h, 3, Dz);

    Eigen::SparseMatrix<double> t;
    igl::cat(1, Dx, Dy, t);
    igl::cat(1, t, Dz, G);
}
