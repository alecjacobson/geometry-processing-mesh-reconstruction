#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

	// Note: DO NO USE D.insert() for sparse matrices its super slow 
	// Use this Triplet thing to the non zero elements in a list, and then 
	// send that to the sparse matrix to read from. MUCH FASTER.
	// TODO: change fd_interpolate to use Triplets

	//Setup the proper dimensions for D, depending on direction
	int sx, sy, sz;
	if (dir == 0)	   { D.resize((nx - 1)*ny*nz, nx*ny*nz); sx = nx - 1; sy = ny; sz = nz; }
	else if (dir == 1) { D.resize(nx*(ny-1)*nz, nx*ny*nz); sx = nx; sy = ny - 1; sz = nz;}
	else if (dir == 2) { D.resize(nx*ny*(nz-1), nx*ny*nz); sx = nx; sy = ny; sz = nz -1;}

	typedef Eigen::Triplet<double> Triple;
	std::vector<Triple> triplets;
	
	//Reserve the two elements of the derivative for each staggered grid location
	triplets.reserve(sx * sy * sz * 2);  

	//Loop through staggered grid and populate derivative values on primary grid.
	for (int x = 0; x < sx; x++) {
		for (int y = 0; y < sy; y++) {
			for (int z = 0; z < sz; z++) {

				//At the grid position i.e l = i-1 -> -1/h
				triplets.push_back(Triple(x + y*sx + z*sx*sy, x + y*nx + z*nx*ny, -1/h));
								
				if (dir == 0) {
					triplets.push_back(Triple(x + y*sx + z*sx*sy, (x + 1) + y*nx + z*nx*ny, 1 / h));
				}
				else if (dir == 1) {
					triplets.push_back(Triple(x + y*sx + z*sx*sy, x + (y+1)*nx + z*nx*ny, 1 / h));
				}
				else if (dir == 2) {
					triplets.push_back(Triple(x + y*sx + z*sx*sy, x + y*nx + (z + 1)*nx*ny, 1 / h));
				}
			}
		}
	}

	D.setFromTriplets(triplets.begin(), triplets.end());

}
