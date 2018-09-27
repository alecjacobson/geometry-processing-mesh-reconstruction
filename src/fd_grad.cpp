#include "fd_grad.h"
#include "igl/cat.h"
#include "fd_partial_derivative.h"
#include <iostream>

using namespace Eigen;

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  SparseMatrix<double> Dx, Dy, Dz, buffer;
  fd_partial_derivative(nx, ny, nz, h, 0, Dx);
  fd_partial_derivative(nx, ny, nz, h, 1, Dy);
  fd_partial_derivative(nx, ny, nz, h, 2, Dz);
  igl::cat(1, Dx, Dy, buffer);
  igl::cat(1, buffer, Dz, G);
}
