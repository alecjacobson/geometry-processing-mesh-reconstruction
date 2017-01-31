#include "fd_grad.h"
#include "fd_partial_derivative.h"

#include <iostream>

/*
01234567801234567801234567801234567801234567801234567801234567801234567
 */


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

  // std::cout << "Dx" << std::endl;
  // std::cout << Dx << std::endl << std::endl;
  // std::cout << "Dy" << std::endl;
  // std::cout << Dy << std::endl << std::endl;
  // std::cout << "Dz" << std::endl;
  // std::cout << Dz << std::endl << std::endl;
  
  G.reserve( Dx.nonZeros() + Dy.nonZeros() + Dz.nonZeros() );
  std::vector< Eigen::Triplet<double> > t;
  extractTriplets( Dx, t, 0 );
  extractTriplets( Dy, t, Dx.rows() );
  extractTriplets( Dz, t, Dx.rows()+Dy.rows() );
  G.setFromTriplets( t.begin(), t.end() );
  G.finalize();

  //std::cout << "D's are size " << Dx.rows() << "x" << Dx.cols() << std::endl;
  //std::cout << "After concatenating, G is " << G.rows() << "x" << G.cols() << std::endl;
}

/*
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


  std::cout << "Dx rows" << Dx.rows() << std::endl;
  std::cout << "Dy rows" << Dy.rows() << std::endl;
  std::cout << "Dz rows" << Dz.rows() << std::endl;
  std::cout << "Dx cols" << Dx.cols() << std::endl;
  std::cout << "Dy cols" << Dy.cols() << std::endl;
  std::cout << "Dz cols" << Dz.cols() << std::endl;
  
  G.reserve( Dx.nonZeros() + Dy.nonZeros() + Dz.nonZeros() );
  std::cout << "Invoking for columns= " << Dx.cols() << std::endl;
  for( unsigned int c=0; c<Dx.cols(); ++c )
  {
      for( Eigen::SparseMatrix<double>::InnerIterator itDx(Dx, c); itDx; ++itDx )
          G.insertBack( itDx.row(), c) = itDx.value();
      for( Eigen::SparseMatrix<double>::InnerIterator itDy(Dy, c); itDy; ++itDy )
          G.insertBack( itDy.row(), c) = itDy.value();
      for( Eigen::SparseMatrix<double>::InnerIterator itDz(Dz, c); itDz; ++itDz )
          G.insertBack( itDz.row(), c) = itDz.value();
  }
  G.finalize();

  std::cout << "D's are size " << Dx.rows() << "x" << Dx.cols() << std::endl;
  std::cout << "After concatenating, G is " << G.rows() << "x" << G.cols() << std::endl;
}
*/
