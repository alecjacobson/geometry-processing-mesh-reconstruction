#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <Eigen/Core>

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> &G) {
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  Eigen::SparseMatrix<double> dx, dy, dz;

  fd_partial_derivative(nx, ny, nz, h, 0, dx);
  fd_partial_derivative(nx, ny, nz, h, 1, dy);
  fd_partial_derivative(nx, ny, nz, h, 2, dz);

  G.resize(dx.rows() + dy.rows() + dz.rows(), dx.cols());
  G.setZero();

  std::vector<Eigen::Triplet<double>> tripletList;
  for (int k = 0; k < dx.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(dx, k); it; ++it) {
      tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
    }
  }
  for (int k = 0; k < dy.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(dy, k); it; ++it) {
      tripletList.push_back(Eigen::Triplet<double>(it.row() + dx.rows(), it.col(), it.value()));
    }
  }
  for (int k = 0; k < dz.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(dz, k); it; ++it) {
      tripletList.push_back(Eigen::Triplet<double>(it.row() + dx.rows() + dy.rows(), it.col(), it.value()));
    }
  }
  G.setFromTriplets(tripletList.begin(), tripletList.end());
  ////////////////////////////////////////////////////////////////////////////
}
