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
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
	W.setZero();
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(8 * P.rows());

	for (int r = 0; r < P.rows(); r++) {
		double x = P(r, 0);
		double y = P(r, 1);
		double z = P(r, 2);

		int i, j, k;
		// x-component
		i = floor((x - corner(0)) / h);
		j = floor((y - corner(1)) / h);
		k = floor((z - corner(2)) / h);

		// The coordinates of the left-bottom-front
		// point of the grid surrounding P
		double grid_corner_x, grid_corner_y, grid_corner_z;
		grid_corner_x = corner(0) + h * i;
		grid_corner_y = corner(1) + h * j;
		grid_corner_z = corner(2) + h * k;

		// The fraction of the coordinate of P relative to the left-bottom-front
		// point of the surrounding grid
		double tx, ty, tz;
		tx = (x - grid_corner_x) / h;
		ty = (y - grid_corner_y) / h;
		tz = (z - grid_corner_z) / h;
		triplets.push_back({r, i + j * nx + k * nx * ny, (1 - tx) * (1 - ty) * (1 - tz)});
		triplets.push_back({r, i + (j + 1) * nx + k * ny * nx, (1 - tx) * ty * (1 - tz)});
		triplets.push_back({r, i + j * nx + (k + 1) * ny * nx, (1 - tx) * (1 - ty) * (tz)});
		triplets.push_back({r, i + (j + 1) * nx + (k + 1) * ny * nx, (1 - tx) * (ty) * (tz)});
		triplets.push_back({r, i + 1 + (j)* nx + (k)* ny * nx, tx * (1 - ty) * (1 - tz)});
		triplets.push_back({r, i + 1 + (j + 1) * nx + (k)* ny * nx, (tx) * (ty) * (1 - tz)});
		triplets.push_back({r, i + 1 + (j)* nx + (k + 1) * ny * nx, (tx) * (1 - ty) * (tz)});
		triplets.push_back({r, i + 1 + (j + 1) * nx + (k + 1) * ny * nx, (tx) * (ty) * (tz)});
	}
	W.setFromTriplets(triplets.begin(), triplets.end());

  ////////////////////////////////////////////////////////////////////////////
}
