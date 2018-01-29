#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  int dx = (nx-1)*ny*nz;
  int dy = nx*(ny-1)*nz;
  int dz = nx*ny*(nz-1);

  Eigen::SparseMatrix<double> Dx = Eigen::SparseMatrix<double>(dx, nx*ny*nz);
  Eigen::SparseMatrix<double> Dy = Eigen::SparseMatrix<double>(dy, nx*ny*nz);
  Eigen::SparseMatrix<double> Dz = Eigen::SparseMatrix<double>(dz, nx*ny*nz);

  fd_partial_derivative(nx, ny, nz, h, 0, Dx);
  fd_partial_derivative(nx, ny, nz, h, 1, Dy);
  fd_partial_derivative(nx, ny, nz, h, 2, Dz);

  // G.block(0, 0, dx, nx*ny*nz) = Dx;
  // G.block(dx, nx*ny*nz, dx+dy, nx*ny*nz) = Dy;
  // G.block(dx+dy, nx*ny*nz, dx+dy+dz, nx*ny*nz) = Dz;
  ////////////////////////////////////////////////////////////////////////////
}
