#include "fd_partial_derivative.h"

typedef Eigen::Triplet<double> T;
extern int gidx(double i, double j, double k, int nx, int ny);
extern void intervalue2(std::vector<T> &coef, int row, int aidx, int bidx,
		int dir);

void fd_partial_derivative(const int nx, const int ny, const int nz,
		const double h, const int dir, Eigen::SparseMatrix<double> & D) {
	////////////////////////////////////////////////////////////////////////////
	// Add your code here
	int nnum = nx * ny * nz;
	D(nnum * 3, nnum * 3);

	// first order
	std::vector<T> coef;
	for (int i = 0; i < nx - 1; i++)
		for (int j = 0; j < ny - 1; j++)
			for (int k = 0; k < nz - 1; k++) {
				int st = gidx(i, j, k, nx, ny);
				int edx = gidx(i + 1, j, k, nx, ny);
				int edy = gidx(i, j + 1, k, nx, ny);
				int edz = gidx(i, j, k + 1, nx, ny);

				int row = st * 1;
				if (dir == 0) {
					intervalue2(coef, row, st, edx, 0);
				} else {
					if (dir == 1)
						intervalue2(coef, row, st, edy, 1);
					else
						intervalue2(coef, row, st, edz, 2);
				}
			}

	D.setFromTriplets(coef.begin(), coef.end());
	////////////////////////////////////////////////////////////////////////////
}
