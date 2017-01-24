#include "fd_interpolate.h"

#include <cmath>
#include <algorithm>

void fd_interpolate(
	const int nx,
	const int ny,
	const int nz,
	const double h,
	const Eigen::RowVector3d & corner,
	const Eigen::MatrixXd & P,
	Eigen::SparseMatrix<double> & W)
{
	W.resize(P.rows(), nx*ny*nz);
	std::vector<Eigen::Triplet<double>> entries;
	entries.reserve(8 * P.rows());

	for (int i = 0; i < P.rows(); i++) {
		auto grid_point = (P.row(i) - corner) / h;

		double weights[3][2];
		Eigen::RowVector3i index;
		for (int j = 0; j < 3; j++) {
			double integral;
			double frac = std::modf(grid_point(j), &integral);

			index(j) = (int)integral;
			weights[j][0] = 1.0 - frac;
			weights[j][1] = frac;
		}

		for (int x = 0; x <= 1; x++) {
			for (int y = 0; y <= 1; y++) {
				for (int z = 0; z <= 1; z++) {
					const int col = index.x() + x + (index.y() + y)*nx + (index.z() + z)*nx*ny;
					entries.emplace_back(
						i,
						col,
						weights[0][x] * weights[1][y] * weights[1][z]);
				}
			}
		}
	}

	W.setFromTriplets(entries.cbegin(), entries.cend());
}

