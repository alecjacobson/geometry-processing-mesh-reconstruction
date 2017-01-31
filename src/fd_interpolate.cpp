#include "fd_interpolate.h"

#include <iostream>

namespace {
double
myFloor( double a )
{
    return (int)floor(a);
}

// Since our grid is given by nx,ny,nz, with cell size h, we have a cubic grid.
// Every cube in our grid will get the same pattern of the contribution of the
// interpolated weights.  The weights themselves for each of the 3 cube dimensions
// are w and 1-w respectively, where w is the normalized distance from the origin
// of the cube. If the cube coordinates are given as x,y,z with lower indices _0
// and _1, ; then the origin or the cube is (lower left hand side)
// origin: (x_0, y_0, z_0)
// upper right back-most point is: (x_1, y_1, z_1).
// So we can store the weights for each axis, for each index as:
// index 0 -> w, and index 1 -> (1-w)
// There are six values, 2 values for each dimension (3d). So for all eight corner
// points we get 2^3 combinations = 8 points of the cube
Eigen::Matrix<double,2,3>
getCubeWeights( Eigen::RowVector3d &pGrid )
{
    Eigen::Matrix<double,2,3> C;
    Eigen::RowVector3d w;

    // from: http://stackoverflow.com/questions/33786662/apply-function-to-all-eigen-matrix-element
    w = pGrid.unaryExpr( &myFloor );
    w = pGrid - w;
    for( int i=0; i<2; ++i )
    {
        for( int j=0; j<3; ++j )
        {
            C(i,j) = ( i ? w(j) : (1-w(j)) );
        }
    }

    //std::cout << "Point " << pGrid << " gets w " << w << " and C is " << C << std::endl;

    return C;
}

// given i,j,k, return the computed index as given by:
// g(i+j*n_x+k*n_y*n_x)
// for the staggered grid subscripts i−½,j,ki−½,j,k we will assume that
// Dxi−½,j,k(ℓ)Di−½,j,kx(ℓ) refers to the matrix entry
// Dx(i+j*n_x+k*n_y*n_x,l), where the i−½i−½ has been rounded down.
int
getIndex( int i, int j, int k, int nx, int ny, int nz )
{
    return( i + nx*( j + k*ny ) );
}

int
convertLocationToIdx( int x, int y, int z,
                      int i, int j, int k,
                      int nx, int ny, int nz )
{
    return getIndex( x+i, y+j, z+k, nx, ny, nz );
}

} // anonymous namespace



void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
    ////////////////////////////////////////////////////////////////////////////
    // Add your code here
    ////////////////////////////////////////////////////////////////////////////

    // construct a sparse matrix W of trilinear interpolation weights so that
    // P = W * x
    std::vector< Eigen::Triplet<double> > t;

    Eigen::Vector3i sz = {nx, ny, nz};
    const int nP = P.rows();
    W.resize( nP, nx*nz*ny );
    for( int l = 0; l < nP; ++l )
    {
        Eigen::RowVector3d p = P.row(l);
        Eigen::RowVector3d pInGrid = (p - corner);
        Eigen::RowVector3d pGrid = pInGrid / h; // normalized grid coords;

        Eigen::RowVector3i gridIdx;
        for( int q = 0; q < 3; q++ )
        {
            gridIdx(q) = (int)std::floor( pGrid(q) );
            // clamp to grid...
            assert( gridIdx(q) < (sz(q)-1) );
            assert( gridIdx(q) >= 0 );
            if( gridIdx(q) > (sz(q)-1) )
                gridIdx(q) = sz(q);
            if( gridIdx(q) < 0 )
                gridIdx(q) = 0;
        }

        Eigen::Matrix<double,2,3> Cw = getCubeWeights( pGrid );

        // do all 8 cube points using our Cw
        for( int i=0; i < 2; ++i )
            for( int j=0; j < 2; ++j )
                for( int k=0; k<2; ++k )
                {
                    int idx = convertLocationToIdx( gridIdx(0), gridIdx(1), gridIdx(2),
                                                    i, j, k,
                                                    nx, ny, nz );
                    t.push_back( Eigen::Triplet<double>{ l, idx,
                                                         Cw( i, 0 )*
                                                         Cw( j, 1 )*
                                                         Cw( k, 2 ) } );
                }
    }

    // Fill in our sparse eigen matrix with the triplets as in the Eigen example...
    W.setFromTriplets( t.begin(), t.end() );
}
