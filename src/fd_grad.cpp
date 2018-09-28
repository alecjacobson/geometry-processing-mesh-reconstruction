#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <igl/cat.h>
#include <iostream>

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  // Construct gradient:
  Eigen::SparseMatrix<double> Dx, Dy, Dz, Inter;
  fd_partial_derivative(nx, ny, nz, h, 0, Dx);
  fd_partial_derivative(nx, ny, nz, h, 1, Dy);
  fd_partial_derivative(nx, ny, nz, h, 2, Dz);

  // Cat:
  Inter.resize(Dx.rows() + Dy.rows(), Dx.cols());
  igl::cat(1, Dx, Dy, Inter);
  G.resize(Inter.rows() + Dz.rows(), Dx.cols());
  igl::cat(1, Inter, Dz, G);
  G.makeCompressed();
}
