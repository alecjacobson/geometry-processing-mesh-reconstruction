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
  Eigen::SparseMatrix<double> D1;
  Eigen::SparseMatrix<double> D2;
  Eigen::SparseMatrix<double> D3;

  fd_partial_derivative(nx, ny, nz, h, 0, D1);
  fd_partial_derivative(nx, ny, nz, h, 1, D2);
  fd_partial_derivative(nx, ny, nz, h, 2, D3);

  G = igl::cat(1, D1, D2);
  G = igl::cat(1, G, D3);
}
