#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <iostream>

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
    int m = nx*ny*nz;
    Eigen::SparseMatrix<double> Dx((nx-1)*ny*nz,m), Dy((ny-1)*nx*nz,m), Dz((nz-1)*ny*nx,m);
    fd_partial_derivative(nx,ny,nz,h,0,Dx);
    fd_partial_derivative(nx,ny,nz,h,1,Dy);
    fd_partial_derivative(nx,ny,nz,h,2,Dz);
    
    Eigen::SparseMatrix<double> C(Dx.rows()+Dy.rows(),Dx.cols());
    igl::cat(1,Dx,Dy,C);
    igl::cat(1,C,Dz,G);
}
