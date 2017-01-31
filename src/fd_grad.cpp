#include "fd_grad.h"
#include "fd_partial_derivative.h"

#include <iostream>

void
extractTriplets( Eigen::SparseMatrix<double> &a,
                 std::vector< Eigen::Triplet<double> > &t,
                 int rowOffset = 0 )
{
  for( unsigned int i = 0; i < a.outerSize(); ++i )
      for( Eigen::SparseMatrix<double>::InnerIterator it(a, i); it; ++it )
          t.push_back( Eigen::Triplet<double>( it.row()+rowOffset, it.col(),
                                               it.value() ) );
}

// Given a regular finite-difference grid described by the number of nodes
// on each side (nx, ny and nz), and the grid spacing (h), construct a
// sparse matrix G to compute gradients with each component on its
// respective staggered grid.

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  Eigen::SparseMatrix<double> Dx, Dy, Dz;
  fd_partial_derivative( nx, ny, nz, h, 0, Dx );
  fd_partial_derivative( nx, ny, nz, h, 1, Dy );
  fd_partial_derivative( nx, ny, nz, h, 2, Dz );

  G.resize( Dx.rows() + Dy.rows() + Dz.rows(),
            Dx.cols() );


  std::cout << "Dx " << Dx.rows() << "x" << Dx.cols() << std::endl;
  std::cout << "Dy " << Dy.rows() << "x" << Dy.cols() << std::endl;
  std::cout << "Dz " << Dz.rows() << "x" << Dz.cols() << std::endl;

  G.reserve( Dx.nonZeros() + Dy.nonZeros() + Dz.nonZeros() );
  std::vector< Eigen::Triplet<double> > t;
  extractTriplets( Dx, t, 0 );
  extractTriplets( Dy, t, Dx.rows() );
  extractTriplets( Dz, t, Dx.rows()+Dy.rows() );

  G.setFromTriplets( t.begin(), t.end() );
  G.finalize();
}
