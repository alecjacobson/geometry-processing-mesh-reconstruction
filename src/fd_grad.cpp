#include "fd_grad.h"
#include "fd_partial_derivative.h"

#include "igl/cat.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  Eigen::SparseMatrix<double> Dx;
  fd_partial_derivative(nx, ny, nz, h, 0, Dx);

  Eigen::SparseMatrix<double> Dy;
  fd_partial_derivative(nx, ny, nz, h, 1, Dy);

  Eigen::SparseMatrix<double> Dz;
  fd_partial_derivative(nx, ny, nz, h, 2, Dz);

  Eigen::SparseMatrix<double> T;
  igl::cat(1, Dx, Dy, T);
  igl::cat(1, T, Dz, G);
}
