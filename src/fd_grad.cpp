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

  //concatenate method copy from stackoverflow
  G.reserve(Dx.nonZeros() + Dy.nonZeros() + Dz.nonZeros());
  for(int c=0; c<Dx.cols(); ++c) {
  for(Eigen::SparseMatrix<double>::InnerIterator itDx(Dx, c); itDx; ++itDx)
    G.insertBack(itDx.row(), c) = itDx.value();
  for(Eigen::SparseMatrix<double>::InnerIterator itDy(Dy, c); itDy; ++itDy)
    G.insertBack(itDy.row(), c) = itDy.value();
  for(Eigen::SparseMatrix<double>::InnerIterator itDz(Dz, c); itDz; ++itDz)
    G.insertBack(itDz.row(), c) = itDz.value();
  }

  G.finalize();

  ////////////////////////////////////////////////////////////////////////////
}
