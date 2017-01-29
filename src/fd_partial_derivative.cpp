#include "fd_partial_derivative.h"
#include <vector>

typedef Eigen::Triplet<double> tri;

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  // Computing the partial derivative in the direction "dir", x = 0; y = 1; z = 2;
  // using the finite difference method, all in the staggered grid
  //                    / -1  if ℓ = i-1
  // Dx_{i-½,j,k}(ℓ) = |   1  if ℓ = i
  //                    \  0  otherwise
  // Dx is of size (nx-1)nynz ✖ nxnynz

  int32_t Nx = nx, Ny = ny, Nz = nz;
  // the offset is the adjacent cell in the matrix depending on direction
  int32_t offset = 0;
  switch(dir)
  {
  case 0:
    Nx--;
    offset = 1;
    break;
  case 1:
    Ny--;
    offset = nx;
    break;
  case 2:
    Nz--;
    offset = nx*ny;
    break;
  }

  // size the matrix according to direction, and has 2 entries per row
  D.resize(Nx*Ny*Nz , nx*ny*nz);
  std::vector<tri> values;
  values.reserve(2 * Nx*Ny*Nz);

  for( int32_t i = 0 ; i < Nx ; i++ )
  {
    for( int32_t j = 0 ; j < Ny ; j++ )
    {
      for( int32_t k = 0 ; k < Nz ; k++ )
      {
	values.push_back(tri(i + Nx*j + Nx*Ny*k, i + nx*j + nx*ny*k,          -1 ));
	values.push_back(tri(i + Nx*j + Nx*Ny*k, i + nx*j + nx*ny*k + offset,  1 ));
      }
    }
  }

  D.setFromTriplets(values.begin(), values.end());
}
