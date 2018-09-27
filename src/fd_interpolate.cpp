#include "fd_interpolate.h"

typedef Eigen::Triplet<double> T;

void intervalue(std::vector<T> &coef, int aidx, int bidx, int value) {
	T tmp;
	for (int i = 0; i < 3; i++) {
		tmp[0] = aidx * 3 + i;
		tmp[1] = bidx * 3 + i;
		tmp[2] = value;
		coef.push_back(tmp);
	}
}

int gidx(double i, double j, double k, int nx, int ny) {
	return i + nx * (j + k * ny);
}

void fd_interpolate(const int nx, const int ny, const int nz, const double h,
		const Eigen::RowVector3d & corner, const Eigen::MatrixXd & P,
		Eigen::SparseMatrix<double> & W) {
	////////////////////////////////////////////////////////////////////////////
	// Add your code here

	int pnum = P.rows();
	int nnum = nx * ny * nz;

	Eigen::MatrixXd P2(pnum * 3, 1);
	for (int i = 0; i < pnum; i++) {
		for (int j = 0; j < 3; j++) {
			P2(i * 3 + j, 0) = P(i, j);
		}
	}

	// size of w should be pnum * 3 x nnum * 3

	std::vector<T> coef;
	for (int p = 0; p < pnum; p++) {
		double px = P(p, 0) - corner[0];
		double py = P(p, 1) - corner[1];
		double pz = P(p, 2) - corner[2];

		double xind = px / h;
		double yind = py / h;
		double zind = pz / h;

		// 8 equations
		double xfloor = std::floor(px);
		double yfloor = std::floor(py);
		double zfloor = std::floor(pz);
		double xceil = xfloor + 1;
		double yceil = yfloor + 1;
		double zceil = zfloor + 1;

		// left, left, left
		double idx = gidx(xfloor, yfloor, zfloor, nx, ny);
		double value = (xceil - px) * (yceil - py) * (zceil - pz);
		intervalue(coef, p, idx, value);

		idx = gidx(xfloor, yceil, zfloor, nx, ny);
		value = (xceil - px) * (py - yfloor) * (zceil - pz);
		intervalue(coef, p, idx, value);

		idx = gidx(xceil, yceil, zfloor, nx, ny);
		value = (px - xfloor) * (py - yfloor) * (zceil - pz);
		intervalue(coef, p, idx, value);

		idx = gidx(xceil, yfloor, zfloor, nx, ny);
		value = (px - xfloor) * (yceil - py) * (zceil - pz);
		intervalue(coef, p, idx, value);

		///////////////////////////////////////////////////////
		idx = gidx(xfloor, yfloor, zceil, nx, ny);
		value = (xceil - px) * (yceil - px) * (pz - zfloor);
		intervalue(coef, p, idx, value);

		idx = gidx(xfloor, yceil, zceil, nx, ny);
		value = (xceil - px) * (px - yfloor) * (pz - zfloor);
		intervalue(coef, p, idx, value);

		idx = gidx(xceil, yceil, zceil, nx, ny);
		value = (px - xfloor) * (px - yfloor) * (pz - zfloor);
		intervalue(coef, p, idx, value);

		idx = gidx(xceil, yfloor, zceil, nx, ny);
		value = (px - xfloor) * (yceil - py) * (pz - zfloor);
		intervalue(coef, p, idx, value);
	}

	W.setFromTriplets(coef.begin(), coef.end());

	////////////////////////////////////////////////////////////////////////////
}
