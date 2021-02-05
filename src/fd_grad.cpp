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
  Eigen::SparseMatrix<double>
  Gx((nx-1)*ny*nz, nx * ny * nz),
  Gy(nx*(ny-1)*nz, nx * ny * nz),
  Gz(nx*ny*(nz-1), nx * ny * nz),
  result;

  fd_partial_derivative(nx, ny, nz, h, 0, Gx);
  fd_partial_derivative(nx, ny, nz, h, 1, Gy);
  fd_partial_derivative(nx, ny, nz, h, 2, Gz);

  igl::cat(1, Gx, Gy, result);
  igl::cat(1, result, Gz, G);
}
