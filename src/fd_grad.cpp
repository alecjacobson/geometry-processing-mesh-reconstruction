#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  int volume = nx*ny*nz;
  Eigen::SparseMatrix<double> Gx((nx-1)*ny*nz, volume);
  Eigen::SparseMatrix<double> Gy(nx*(ny-1)*nz, volume);
  Eigen::SparseMatrix<double> Gz(nx*ny*(nz-1), volume);

  fd_partial_derivative(nx, ny, nz, h, 0, Gx);
  fd_partial_derivative(nx, ny, nz, h, 1, Gy);
  fd_partial_derivative(nx, ny, nz, h, 2, Gz);

  Eigen::SparseMatrix <double> tmp;

  // Why is dim=1 or 2 and not 0 or 1 in igl::cat????
  igl::cat(1, Gx, Gy, tmp);
  igl::cat(1, tmp, Gz, G);
}
