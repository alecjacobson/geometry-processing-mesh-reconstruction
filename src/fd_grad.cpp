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
  Eigen::SparseMatrix<double> G1,G2,G3;
  fd_partial_derivative(nx,ny,nz,h,0,G1);
  fd_partial_derivative(nx,ny,nz,h,1,G2);
  fd_partial_derivative(nx,ny,nz,h,2,G3);
  G=igl::cat(1,G1,G2);
  G=igl::cat(1,G,G3);
}
