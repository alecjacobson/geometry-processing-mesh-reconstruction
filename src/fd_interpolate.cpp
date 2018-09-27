#include "fd_interpolate.h"
#include <assert.h>

typedef Eigen::Triplet<double> T;

void intervalue(std::vector<T> &coef, double aidx, double bidx, double value) {

	for (int i = 0; i < 3; i++) {
		T tmp(aidx * 3 + i, bidx * 3 + i, value);
		coef.push_back(tmp);
	}
}

double gidx(double i, double j, double k, int nx, int ny, int nz) {
	double idx = i + nx * (j + k * ny);
	assert(idx < (nx - 1) * (ny - 1) * (nz - 1));
	return idx;
}

void fd_interpolate(const int nx, const int ny, const int nz, const double h,
		const Eigen::RowVector3d & corner, const Eigen::MatrixXd & P,
		Eigen::SparseMatrix<double> & W) {
	////////////////////////////////////////////////////////////////////////////
	// Add your code here

	int pnum = P.rows();

	// size of w should be pnum * 3 x (nx-1)(ny-1)(nz-1) * 3
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
		double idx = gidx(xfloor, yfloor, zfloor, nx, ny, nz);
		double value = (xceil - px) * (yceil - py) * (zceil - pz);
		intervalue(coef, p, idx, value);

		idx = gidx(xfloor, yceil, zfloor, nx, ny, nz);
		value = (xceil - px) * (py - yfloor) * (zceil - pz);
		intervalue(coef, p, idx, value);

		idx = gidx(xceil, yceil, zfloor, nx, ny, nz);
		value = (px - xfloor) * (py - yfloor) * (zceil - pz);
		intervalue(coef, p, idx, value);

		idx = gidx(xceil, yfloor, zfloor, nx, ny, nz);
		value = (px - xfloor) * (yceil - py) * (zceil - pz);
		intervalue(coef, p, idx, value);

		///////////////////////////////////////////////////////
		idx = gidx(xfloor, yfloor, zceil, nx, ny, nz);
		value = (xceil - px) * (yceil - px) * (pz - zfloor);
		intervalue(coef, p, idx, value);

		idx = gidx(xfloor, yceil, zceil, nx, ny, nz);
		value = (xceil - px) * (px - yfloor) * (pz - zfloor);
		intervalue(coef, p, idx, value);

		idx = gidx(xceil, yceil, zceil, nx, ny, nz);
		value = (px - xfloor) * (px - yfloor) * (pz - zfloor);
		intervalue(coef, p, idx, value);

		idx = gidx(xceil, yfloor, zceil, nx, ny, nz);
		value = (px - xfloor) * (yceil - py) * (pz - zfloor);
		intervalue(coef, p, idx, value);
	}

	W.setFromTriplets(coef.begin(), coef.end());

	////////////////////////////////////////////////////////////////////////////
}

