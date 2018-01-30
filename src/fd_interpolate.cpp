#include "fd_interpolate.h"

void fd_interpolate(
		const int nx,
		const int ny,
		const int nz,
		const double h,
		const Eigen::RowVector3d & corner,
		const Eigen::MatrixXd & P,
		Eigen::SparseMatrix<double> & W) {
	W.resize(P.rows(), nx * ny * nz);

	auto toGridIndex = [nx, ny, nz](int x, int y, int z) -> int {
		return x + nx * y + nx * ny * z;
	};

	// As the Eigen guide says to, we build the matrix from triplets.
	typedef Eigen::Triplet<double> Triplet;
	std::vector<Triplet> tripletList;

	int rowCount = P.rows();
	for (int r = 0; r < rowCount; r++) {

		// How many h's from the corner the point is.
		double xFromStart = (P(r, 0) - corner(0)) / h;
		double yFromStart = (P(r, 1) - corner(1)) / h;
		double zFromStart = (P(r, 2) - corner(2)) / h;

		// Indices of the lower corner of the grid cell.
		int xGridIndex = floor(xFromStart);
		int yGridIndex = floor(yFromStart);
		int zGridIndex = floor(zFromStart);

		// Where within the current grid cell the point is.
		double x = xFromStart - xGridIndex;
		double y = yFromStart - yGridIndex;
		double z = zFromStart - zGridIndex;

		// The weights computed by hand
		double w0 = (1 - x) * (1 - z) * (1 - y);
		double w1 = (1 - x) * z * (1 - y);
		double w2 = (1 - x) * (1 - z) * y;
		double w3 = (1 - x) * z * y;
		double w4 = x * (1 - z) * (1 - y);
		double w5 = x * z * (1 - y);
		double w6 = x * (1 - z) * y;
		double w7 = x * z * y;

		tripletList.push_back(Triplet(r, toGridIndex(xGridIndex, yGridIndex, zGridIndex), w0));
		tripletList.push_back(Triplet(r, toGridIndex(xGridIndex, yGridIndex, zGridIndex + 1), w1));
		tripletList.push_back(Triplet(r, toGridIndex(xGridIndex, yGridIndex + 1, zGridIndex), w2));
		tripletList.push_back(Triplet(r, toGridIndex(xGridIndex, yGridIndex + 1, zGridIndex + 1), w3));

		tripletList.push_back(Triplet(r, toGridIndex(xGridIndex + 1, yGridIndex, zGridIndex), w4));
		tripletList.push_back(Triplet(r, toGridIndex(xGridIndex + 1, yGridIndex, zGridIndex + 1), w5));
		tripletList.push_back(Triplet(r, toGridIndex(xGridIndex + 1, yGridIndex + 1, zGridIndex), w6));
		tripletList.push_back(Triplet(r, toGridIndex(xGridIndex + 1, yGridIndex + 1, zGridIndex + 1), w7));
	}

	W.setFromTriplets(tripletList.begin(), tripletList.end());
}
