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
  int dx = (nx-1)*ny*nz;
  int dy = nx*(ny-1)*nz;
  int dz = nx*ny*(nz-1);

  Eigen::SparseMatrix<double> Dx = Eigen::SparseMatrix<double>(dx, nx*ny*nz);
  Eigen::SparseMatrix<double> Dy = Eigen::SparseMatrix<double>(dy, nx*ny*nz);
  Eigen::SparseMatrix<double> Dz = Eigen::SparseMatrix<double>(dz, nx*ny*nz);

  fd_partial_derivative(nx, ny, nz, h, 0, Dx);
  fd_partial_derivative(nx, ny, nz, h, 1, Dy);
  fd_partial_derivative(nx, ny, nz, h, 2, Dz);
  stack_sparse_matric(Dx, Dy, Dz, G);
  ////////////////////////////////////////////////////////////////////////////
}

//concatenate method copy from stackoverflow
void stack_sparse_matric(
  Eigen::SparseMatrix<double> &A,
  Eigen::SparseMatrix<double> &B,
  Eigen::SparseMatrix<double> &C, 
  Eigen::SparseMatrix<double> &out) 
{
  std::vector<Eigen::Triplet<double> > tripletList;
  for (int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
      tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
    }
  }
  for (int k = 0; k < B.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
      tripletList.push_back(Eigen::Triplet<double>(A.rows() + it.row(), it.col(), it.value()));
    }
  }
  for (int k = 0; k < C.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
      tripletList.push_back(Eigen::Triplet<double>(A.rows() + B.rows() + it.row(), it.col(), it.value()));
    }
  }
  out.setFromTriplets(tripletList.begin(), tripletList.end());
}
