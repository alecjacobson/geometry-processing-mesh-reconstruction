#include "fd_grad.h"

typedef Eigen::Triplet<double> T;

void intervalue2(std::vector<T> &coef, int row, int aidx, int bidx, int dir) {
	T tmp;
	tmp[0] = row;
	tmp[1] = bidx * 3 + dir;
	tmp[2] = 1;
	coef.push_back(tmp);

	tmp[0] = row;
	tmp[1] = aidx * 3 + dir;
	tmp[2] = -1;
	coef.push_back(tmp);
}

extern int gidx(double i, double j, double k, int nx, int ny);

void fd_grad(const int nx, const int ny, const int nz, const double h,
		Eigen::SparseMatrix<double> & G) {
	////////////////////////////////////////////////////////////////////////////
	// Add your code here
	int nnum = nx * ny * nz;
	G(nnum * 3, nnum * 3);

	// first order
	std::vector<T> coef;
	for (int i = 0; i < nx - 1; i++)
		for (int j = 0; j < ny - 1; j++)
			for (int k = 0; k < nz - 1; k++) {
				int st = gidx(i, j, k, nx, ny);
				int edx = gidx(i + 1, j, k, nx, ny);
				int edy = gidx(i, j + 1, k, nx, ny);
				int edz = gidx(i, j, k + 1, nx, ny);

				int row = st * 3;
				intervalue2(coef, row, st, edx, 0);
				intervalue2(coef, row, st+1, edy, 1);
				intervalue2(coef, row, st+2, edz, 2);
			}
	G.setFromTriplets(coef.begin(), coef.end());
	////////////////////////////////////////////////////////////////////////////
}
