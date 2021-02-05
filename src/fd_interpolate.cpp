#include "fd_interpolate.h"
#include <math.h>

void fd_interpolate(
	const int nx,
	const int ny,
	const int nz,
	const double h,
	const Eigen::RowVector3d & corner,
	const Eigen::MatrixXd & P,
	Eigen::SparseMatrix<double> & W)
{
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(P.rows() * 8);
	for(int i = 0; i < P.rows(); i++) {
		// relative position to corner
		double x = (P(i, 0) - corner[0])/h;
		double y = (P(i, 1) - corner[1])/h;
		double z = (P(i, 2) - corner[2])/h;

		// previous on grid position
		int x1 = floor(x);
		int y1 = floor(y);
		int z1 = floor(z);

		// next on grid position
		int x2 = x1 + 1;
		int y2 = y1 + 1;
		int z2 = z1 + 1;

		// calculate weight for each of eight surrounding point
		tripletList.push_back(T(i, (x1 + y1 * nx + z1 * ny * nx), (x2 - x) * (y2 - y) * (z2 - z)));
		tripletList.push_back(T(i, (x1 + y1 * nx + z2 * ny * nx), (x2 - x) * (y2 - y) * (z - z1)));
		tripletList.push_back(T(i, (x1 + y2 * nx + z1 * ny * nx), (x2 - x) * (y - y1) * (z2 - z)));
		tripletList.push_back(T(i, (x1 + y2 * nx + z2 * ny * nx), (x2 - x) * (y - y1) * (z - z1)));
		tripletList.push_back(T(i, (x2 + y1 * nx + z1 * ny * nx), (x - x1) * (y2 - y) * (z2 - z)));
		tripletList.push_back(T(i, (x2 + y1 * nx + z2 * ny * nx), (x - x1) * (y2 - y) * (z - z1)));
		tripletList.push_back(T(i, (x2 + y2 * nx + z1 * ny * nx), (x - x1) * (y - y1) * (z2 - z)));
		tripletList.push_back(T(i, (x2 + y2 * nx + z2 * ny * nx), (x - x1) * (y - y1) * (z - z1)));
	}
	W.resize(P.rows(), nx * ny * nz);
	W.setFromTriplets(tripletList.begin(), tripletList.end());
}