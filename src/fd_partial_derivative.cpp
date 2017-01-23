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
	D.resize(m, nx*ny*nz);

	std::vector<Eigen::Triplet<double>> entries;
	entries.reserve(2 * m);

	for (int i = 0; i < m; i++) {
		entries.emplace_back(i, i + 1, -1);
		entries.emplace_back(i, i, 1);
	}

	D.setFromTriplets(entries.cbegin(), entries.cend());
}
