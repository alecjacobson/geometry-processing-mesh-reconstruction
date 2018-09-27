#include "fd_partial_derivative.h"

typedef Eigen::Triplet<double> T;
extern int gidx(double i, double j, double k, int nx, int ny, int nz);
extern void intervalue2(std::vector<T> &coef, int row, int aidx, int bidx);

void fd_partial_derivative(const int nx, const int ny, const int nz,
		const double h, const int dir, Eigen::SparseMatrix<double> & D) {
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

				int row = st;
				switch (dir) {
				case 0:
					intervalue2(coef, row, st, edx);
					break;
				case 1:
					intervalue2(coef, row + 1, st, edy);
					break;
				case 2:
					intervalue2(coef, row + 2, st, edz);
					break;
				}
			}

	D.setFromTriplets(coef.begin(), coef.end());
	////////////////////////////////////////////////////////////////////////////
}
