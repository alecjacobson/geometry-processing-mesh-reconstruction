#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  Eigen::SparseMatrix<double> Dx((nx - 1)*ny*nz, nx*ny*nz), Dy(nx*(ny-1)*nz, nx*ny*nz), Dz(nx*ny*(nz-1), nx*ny*nz);
  fd_partial_derivative(nx-1, ny, nz, h, 0, Dx);
  fd_partial_derivative(nx, ny-1, nz, h, 1, Dy);
  fd_partial_derivative(nx, ny, nz-1, h, 2, Dz);
  
  G.reserve(Dx.nonZeros() + Dy.nonZeros() + Dz.nonZeros());
  std::vector<Eigen::Triplet<double> > tripletList;
  for (int k = 0; k < Dx.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Dx, k); it; ++it) {
      tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
    }
  }
  for (int k = 0; k < Dy.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Dy, k); it; ++it) {
      tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
    }
  }
  for (int k = 0; k < Dz.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Dz, k); it; ++it) {
      tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
    }
  }
  G.setFromTriplets(tripletList.begin(), tripletList.end());
}
