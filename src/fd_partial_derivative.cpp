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
	const int bounds[3]{ nx - (dir == 0), ny - (dir == 1), nz - (dir == 2) };
	//Do we need to worry about overflow?
	unsigned int m = bounds[0] * bounds[1] * bounds[2];


	std::vector<Eigen::Triplet<double>> entries;
	entries.reserve(m * 2);

	auto row_index = [&](int x, int y, int z) {
		return x + y*bounds[0] + z*bounds[0] * bounds[1];
	};

	auto col_index = [&](const Eigen::Vector3i& i) {
		return i.x() + i.y()*nx + i.z()*nx*ny;
	};


	for (int x = 0; x < bounds[0]; x++) {
		for (int y = 0; y < bounds[1]; y++) {
			for (int z = 0; z < bounds[2]; z++) {
				Eigen::Vector3i indices(x, y, z);
				entries.emplace_back(row_index(x, y, z), col_index(indices), -1.0);
				indices(dir) += 1;
				entries.emplace_back(row_index(x, y, z), col_index(indices), 1.0);
			}
		}
	}

	D.resize(m, nx*ny*nz);
	D.setFromTriplets(entries.cbegin(), entries.cend());
}
