#include <igl/cat.h>
#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
    Eigen::SparseMatrix<double> Dx;
    Eigen::SparseMatrix<double> Dy;
    Eigen::SparseMatrix<double> Dz;

    fd_partial_derivative(nx, ny, nz, h, X, Dx);
    fd_partial_derivative(nx, ny, nz, h, Y, Dy);
    fd_partial_derivative(nx, ny, nz, h, Z, Dz);

    Eigen::SparseMatrix<double> Dxy;
    igl::cat(1, Dx, Dy, Dxy);
    igl::cat(1, Dxy, Dz, G);
}
