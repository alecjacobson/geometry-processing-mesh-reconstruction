#include "fd_partial_derivative.h"

void fd_partial_derivative(
		const int nx,
		const int ny,
		const int nz,
		const double h,
		const int dir,
		Eigen::SparseMatrix<double> & D) {
	
	// As the Eigen guide says to, we build the matrix from triplets.
	typedef Eigen::Triplet<double> Triplet;
	std::vector<Triplet> tripletList;

	// Iterate over all cells. The limit depends on which direction we're
	// going, so we can avoid going off the end.
	int xLimit = dir == 0 ? nx - 1 : nx;
	int yLimit = dir == 1 ? ny - 1 : ny;
	int zLimit = dir == 2 ? nz - 1 : nz;
	int r = 0;
	for (int x = 0; x < xLimit; x++) {
		for (int y = 0; y < yLimit; y++) {
			for (int z = 0; z < zLimit; z++) {

				// Co-efficient corresponding to one of the points in the
				// finite difference. I divide by h to represent the division
				// by dx.
				int firstIndex = x + y * nx + z * nx * ny;
				tripletList.push_back(Triplet(r, firstIndex, -1.0 / h));

				// Second co-efficient. Again we have to divide by h.
				int secondX = dir == 0 ? x + 1 : x;
				int secondY = dir == 1 ? y + 1 : y;
				int secondZ = dir == 2 ? z + 1 : z;
				int secondIndex = secondX + secondY * nx + 
						secondZ * nx * ny;
				tripletList.push_back(Triplet(r, secondIndex, 1.0 / h));

				r++;
			}
		}
	}

	D.setFromTriplets(tripletList.begin(), tripletList.end());
}
