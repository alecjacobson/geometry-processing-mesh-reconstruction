#include "fd_interpolate.h"

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
	using namespace std;
	typedef Eigen::Triplet<double> T;
	vector<T> triplets;
	triplets.reserve(P.rows()*8*2);
	for (int row = 0; row < P.rows(); row++) {
		double dx = P(row, 0) - corner(0);
		double dy = P(row, 1) - corner(1);
		double dz = P(row, 2) - corner(2);
		int i = int(floor(dx / h));
		double x_value = fmod(dx, h) / h;
		int j = int(floor(dy / h));
		double y_value = fmod(dy, h) / h;
		int k = int(floor(dz / h));
		double z_value = fmod(dz, h) / h;

		// compute the value for the 8 neighboring nodes
		triplets.push_back(T(row, i + nx * (j + k * ny), 
			(1.0 - x_value) * (1.0 - y_value) * (1.0 - z_value)));
		triplets.push_back(T(row, (i + 1) + nx * (j + k * ny),
			x_value * (1.0 - y_value) * (1.0 - z_value)));
		triplets.push_back(T(row, i + nx * ((j + 1) + k * ny),
			(1.0 - x_value) * y_value * (1.0 - z_value)));
		triplets.push_back(T(row, (i + 1) + nx * ((j + 1) + k * ny),
			x_value * y_value * (1.0 - z_value)));
		triplets.push_back(T(row, i + nx * (j + (k + 1) * ny),
			(1.0 - x_value) * (1.0 - y_value) * z_value));
		triplets.push_back(T(row, (i + 1) + nx * (j + (k + 1) * ny),
			x_value * (1.0 - y_value) * z_value));
		triplets.push_back(T(row, i + nx * ((j + 1) + (k + 1) * ny),
			(1.0 - x_value) * y_value * z_value));
		triplets.push_back(T(row, (i + 1) + nx * ((j + 1) + (k + 1) * ny),
			x_value * y_value * z_value));
	}
	W.resize(P.rows(), nx*ny*nz);
	W.setFromTriplets(triplets.begin(), triplets.end());
}
