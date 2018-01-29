#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
	Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(nx*ny*nz*4);
	int x_num = nx, y_num = ny, z_num = nz;
	int ix = 0, jy = 0, kz = 0;
	if (dir == 0) {// partial derivative with regard to x
		x_num--; // set number of grid nodes in direction x
		ix = 1; // set the difference in x direction
	}
	else if (dir == 1) {// partial derivative with regard to y
		y_num--; // set number of grid nodes in direction y
		jy = 1; // set the difference in y direction
	}
	else {// dir == 2 // partial derivative with regard to z
		z_num--; // set number of grid nodes in direction z
		kz = 1; // set the difference in z direction
	}
	double dif = 1.0 / h;
	// compute the partial derivative matrix
	for (int i = 0; i < x_num; i++) {
		for (int j = 0; j < y_num; j++) {
			for (int k = 0; k < z_num; k++) {
				int row = i + x_num * (j + k * y_num);
				int col1 = i + nx * (j + k * ny);
				int col2 = i + ix + nx * ((j + jy) + (k + kz) * ny);
				triplets.push_back(T(row, col1, -dif));
				triplets.push_back(T(row, col2, dif));
			}
		}
	}
	D.resize(x_num*y_num*z_num, nx*ny*nz);
	D.setFromTriplets(triplets.begin(), triplets.end());
}
