#include "fd_grad.h"
#include "fd_partial_derivative.h"

typedef Eigen::Triplet<double> tri;

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  // construct the full gradient matrix from partial directions:
  //      / Dx \
  // G = |  Dy  |
  //      \ Dz /

  Eigen::SparseMatrix<double> Dx, Dy, Dz;

  fd_partial_derivative(nx, ny, nz, h, 0, Dx);
  fd_partial_derivative(nx, ny, nz, h, 1, Dy);
  fd_partial_derivative(nx, ny, nz, h, 2, Dz);

  // now to concat the matrices, sparse matrices do not support block operations :o
  G.resize(Dx.rows() + Dy.rows() + Dz.rows(), Dx.cols());

  std::vector<tri> entries;
  entries.reserve(2 * (Dx.rows() + Dy.rows() + Dz.rows()));

  for( int32_t i = 0 ; i < Dx.outerSize() ; i++ )
    for( Eigen::SparseMatrix<double>::InnerIterator it(Dx, i) ; it ; ++it)
      entries.push_back(tri(it.row(), it.col(), it.value()));

  for( int32_t i = 0 ; i < Dy.outerSize() ; i++ )
    for( Eigen::SparseMatrix<double>::InnerIterator it(Dy, i) ; it ; ++it)
      entries.push_back(tri(Dx.rows() + it.row(), it.col(), it.value()));

  for( int32_t i = 0 ; i < Dz.outerSize() ; i++ )
    for( Eigen::SparseMatrix<double>::InnerIterator it(Dz, i) ; it ; ++it)
      entries.push_back(tri(Dx.rows() + Dy.rows() + it.row(), it.col(), it.value()));

  G.setFromTriplets(entries.begin(), entries.end());
}
