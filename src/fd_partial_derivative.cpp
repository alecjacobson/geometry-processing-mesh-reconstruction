#include "fd_partial_derivative.h"

namespace {
int
getIndex( int i, int j, int k, int nx, int ny, int nz )
{
    return( i + nx*( j + k*ny ) );
}
}

// Given a regular finite-difference grid described by the number of
// nodes on each side (nx, ny and nz), the grid spacing (h), and a
// desired direction, construct a sparse matrix D to compute first
// partial derivatives in the given direction onto the staggered grid
// in that direction.
//
void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
    // D

    // have two grids, the primary and the staggered grid
    // size of the staggered grid is determined from the dir
    int sgnx( nx ), sgny( ny ), sgnz( nz );
    switch( dir )
    {
        case 0:
            sgnx = nx - 1;
            break;
        case 1:
            sgny = ny - 1;
            break;
        case 2:
            sgnz = nz - 1;
            break;
    }


    // Partial is (current - previous)/h
    // division is expensive...so do it once...
    double invH = 1.0 / h;
    
    // D is sparse
    std::vector< Eigen::Triplet<double> > t;

    // D contains 
    // l-th entry in row is -1 for previous value, 1 for current,
    // 0 otherwise
    for( int i = 0; i < sgnx; ++i )
        for( int j = 0; j < sgny; ++j )
            for( int k = 0; k < sgnz; ++k )
            {
                int rowIdx = getIndex( i, j, k, sgnx, sgny, sgnz );

                // previous node
                int colIdx = getIndex( i, j, k, nx, ny, nz );
                t.push_back( Eigen::Triplet<double>{ rowIdx, colIdx,
                                                     -invH } );

                // current node
                switch( dir )
                {
                case 0:
                    colIdx = getIndex( i+1, j, k, nx, ny, nz );
                    t.push_back( Eigen::Triplet<double>{ rowIdx,
                                                         colIdx,
                                                         invH } );
                    break;
                case 1:
                    colIdx = getIndex( i, j+1, k, nx, ny, nz );
                    t.push_back( Eigen::Triplet<double>{ rowIdx,
                                                         colIdx,
                                                         invH } );

                    break;
                case 2:
                    colIdx = getIndex( i, j, k+1, nx, ny, nz );
                    t.push_back( Eigen::Triplet<double>{ rowIdx,
                                                         colIdx,
                                                         invH } );
                    break;
                }
                
            }

    D.resize( sgnx*sgny*sgnz, nx*ny*nz ); // thank god it's sparse...
    D.setFromTriplets( t.begin(), t.end() );
}
