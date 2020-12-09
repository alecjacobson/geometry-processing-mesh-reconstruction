#include "fd_grad.h"

typedef Eigen::Triplet<double> T;

void intervalue2(std::vector<T> &coef, int row, int aidx, int bidx) {
	T tmp(row, aidx, -1);
	T tmp2(row, bidx, 1);
	coef.push_back(tmp);
	coef.push_back(tmp2);
}

extern int gidx(double i, double j, double k, int nx, int ny, int nz);

void fd_grad(const int nx, const int ny, const int nz, const double h,
		Eigen::SparseMatrix<double> & G) {
	////////////////////////////////////////////////////////////////////////////
	// Add your code here

	// first order
	std::vector<T> coef;
	for (int i = 0; i < nx - 1; i++)
		for (int j = 0; j < ny - 1; j++)
			for (int k = 0; k < nz - 1; k++) {

				int st = gidx(i, j, k, nx, ny, nz);
				int edx = gidx(i + 1, j, k, nx, ny, nz);
				int edy = gidx(i, j + 1, k, nx, ny, nz);
				int edz = gidx(i, j, k + 1, nx, ny, nz);

				int row = st * 3;
				intervalue2(coef, row, st, edx);
				intervalue2(coef, row + 1, st, edy);
				intervalue2(coef, row + 2, st, edz);
			}

	G.setFromTriplets(coef.begin(), coef.end());
	////////////////////////////////////////////////////////////////////////////
}

