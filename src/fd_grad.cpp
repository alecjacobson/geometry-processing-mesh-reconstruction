#include "fd_grad.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  int num_cols = nx*ny*nz;
  Eigen::SparseMatrix<double> Dx((nx-1)*ny*nz, num_cols);
  fd_partial_derivative((nx-1),ny,nz,h,0,Dx);
  Eigen::SparseMatrix<double> Dy(nx*(ny-1)*nz, num_cols);
  fd_partial_derivative(nx,(ny-1),nz,h,1,Dy);
  Eigen::SparseMatrix<double> Dz(nx*ny*(nz-1), num_cols);
  fd_partial_derivative(nx,ny,(nz-1),h,2,Dz);

  //G = [Dx;Dy;Dz];
  Eigen::SparseMatrix<double> Dx_Dy((nx-1)*ny*nz + nx*(ny-1)*nz, num_cols);
  igl::cat(1, Dx, Dy, Dx_Dy);
  igl::cat(1, Dx_Dy, Dz, G);
}
