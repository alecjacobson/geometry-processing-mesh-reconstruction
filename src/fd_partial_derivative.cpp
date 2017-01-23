#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
	const int ns[3]{ nx, ny, nz };

	//Do we need to worry about overflow?
	unsigned int m = nx*ny*nz / ns[dir]; //This is equivalent to (nx-1)*ny*nz when dir=0, etc.
	

	std::vector<Eigen::Triplet<double>> entries;
	entries.reserve(m * 2);

	const int bounds[3]{ nx - (dir == 0), ny - (dir == 1), nz - (dir == 2) };

	for (int x = 0; x < bounds[0]; x++) {
		for (int y = 0; y < bounds[1]; y++) {
			for (int z = 0; z < bounds[2]; z++) {
				const int ls[3]{ x,y,z };
				auto index = x + y*bounds[0] + z*bounds[0] * bounds[1];
				entries.emplace_back(index, ls[dir] - 1, -1.0);
				entries.emplace_back(index, ls[dir], -1.0);
			}
		}
	}

	D.resize(m, nx*ny*nz);
	D.setFromTriplets(entries.cbegin(), entries.cend());
}
