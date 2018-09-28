#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include "igl/cat.h"

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

  int m = (nx-1)*ny*nz+ nx*(ny-1)*nz+ nx*ny*(nz-1);
  int n = nx*ny*nz;
  G.resize(m,n);

  Eigen::SparseMatrix<double> Dx,Dy,Dz;
  fd_partial_derivative(nx,ny,nz,h,0,Dx);
  fd_partial_derivative(nx,ny,nz,h,1,Dy);
  fd_partial_derivative(nx,ny,nz,h,2,Dz);
  
  //Then concatenate Dx,Dy,Dz by libigl
  Eigen::SparseMatrix<double> T;
  igl::cat(1,Dx,Dy,T);
  igl::cat(1,T,Dz,G);


}
