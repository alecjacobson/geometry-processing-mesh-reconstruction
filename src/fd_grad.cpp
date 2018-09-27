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
  Eigen::SparseMatrix<double> xmat, ymat, zmat;

  fd_partial_derivative(nx, ny, nz, h, 0, xmat);
  fd_partial_derivative(nx, ny, nz, h, 1, ymat);
  fd_partial_derivative(nx, ny, nz, h, 2, zmat);

  G.resize(xmat.rows() + ymat.rows() + zmat.rows(), xmat.cols());

  G.reserve(xmat.nonZeros() + ymat.nonZeros() + zmat.nonZeros());
  for (int c = 0; c < xmat.cols(); ++c) {
    G.startVec(c);
    for (Eigen::SparseMatrix<double>::InnerIterator itX(xmat, c); itX; ++itX)
      G.insertBack(itX.row(), c) = itX.value();
    for (Eigen::SparseMatrix<double>::InnerIterator itY(xmat, c); itY; ++itY)
      G.insertBack(itY.row() + xmat.rows(), c) = itY.value();
    for (Eigen::SparseMatrix<double>::InnerIterator itZ(xmat, c); itZ; ++itZ)
      G.insertBack(itZ.row() + xmat.rows() + ymat.rows(), c) = itZ.value();
  }
  G.finalize();
  ////////////////////////////////////////////////////////////////////////////
}
