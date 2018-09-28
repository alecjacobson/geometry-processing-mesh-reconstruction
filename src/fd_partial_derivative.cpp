#include "fd_partial_derivative.h"

void fd_partial_derivative(
	const int nx,
	const int ny,
	const int nz,
	const double h,
	const int dir,
	Eigen::SparseMatrix<double> & D)
{
	// compute number of entries
	int x = nx;
	int y = ny;
	int z = nz;
	switch (dir){
		case 0:
			x--;
			break;
		case 1:
			y--;
			break;
		case 2:
			z--;
			break;
	}

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(x * y * z * 2);

	for(int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			for (int k = 0; k < z; k++) {
				int s = i + j * x + k * y * x;
				int t = i + j * nx + k * ny * nx;
				tripletList.push_back(T(s, t, -1/h));
				switch (dir){
					case 0:
						t = (i + 1) + j * nx + k * ny * nx;
						break;
					case 1:
						t = i + (j + 1) * nx + k * ny * nx;
						break;
					case 2:
						t = i + j * nx + (k + 1) * ny * nx;
						break;
				}
				tripletList.push_back(T(s, t, 1/h));
			}
		}
	}
	D.resize(x * y * z, nx * ny * nz);
	D.setFromTriplets(tripletList.begin(), tripletList.end());
}
