#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
    Eigen::SparseMatrix<double> Dx((nx-1)*ny*nz, nx*ny*nz);
    Eigen::SparseMatrix<double> Dy(nx*(ny-1)*nz, nx*ny*nz);
    Eigen::SparseMatrix<double> Dz(nx*ny*(nz-1), nx*ny*nz);

    fd_partial_derivative(nx, ny, nz, h, 0, Dx);
    fd_partial_derivative(nx, ny, nz, h, 1, Dy);
    fd_partial_derivative(nx, ny, nz, h, 2, Dz);

    vstack3(G, Dx, Dy, Dz);
}
