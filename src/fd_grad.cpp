#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <vector>

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  int rows = (nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1);
  int cols = nx*ny*nz;
  G.resize(rows, cols);
  
  Eigen::SparseMatrix<double> Dx, Dy, Dz;
  fd_partial_derivative(nx, ny, nz, h, 0, Dx);
  fd_partial_derivative(nx, ny, nz, h, 1, Dy);
  fd_partial_derivative(nx, ny, nz, h, 2, Dz);
  
  // ideally, there would be a helper for concatenating sparse matrices
  
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  
  for (int i = 0; i < Dx.outerSize(); i++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Dx,i); it; ++it) {
      tripletList.push_back(T(it.row(), it.col(), it.value()));
    }
  }
  for (int i = 0; i < Dy.outerSize(); i++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Dy,i); it; ++it) {
      tripletList.push_back(T(Dx.rows() + it.row(), it.col(), it.value()));
    }
  }
  for (int i = 0; i < Dz.outerSize(); i++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Dz,i); it; ++it) {
      tripletList.push_back(T(Dx.rows() + Dy.rows() + it.row(), it.col(), it.value()));
    }
  }
   
  G.setFromTriplets(tripletList.begin(), tripletList.end());
}
