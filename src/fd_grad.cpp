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
  Eigen::SparseMatrix<double> Dx, Dy, Dz, temp;

  fd_partial_derivative(nx,ny,nz,h,0,Dx);
  fd_partial_derivative(nx,ny,nz,h,1,Dy);
  fd_partial_derivative(nx,ny,nz,h,2,Dz);

  G.resize(Dx.rows() + Dy.rows() + Dz.rows(), nx*ny*nz);

  igl::cat(1, Dx, Dy, temp);
  igl::cat(1, temp, Dz, G);
}