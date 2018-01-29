#include "fd_interpolate.h"
#include <igl/floor.h>
#include <iostream>
#include <typeinfo>

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

	W.resize(P.rows(), nx*ny*nz);

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(P.rows());

	for (int i = 0; i < P.rows(); i++)
	{
		Eigen::RowVector3d point(P(i, 0), P(i,1), P(i,2));
		Eigen::RowVector3d grid_index = (point - corner) / h;
		Eigen::RowVector3d floored_index;
		igl::floor(grid_index, floored_index);

		// the three weights are how much "over" the integer grid values the double values are AKA the decimal part of those numbers
		Eigen::RowVector3d weights = ((grid_index + corner) - floored_index);

		// now we construct all eight
		int x = floored_index(0);
		int y = floored_index(1);
		int z = floored_index(2);

		double w_x = weights(0);
		double w_y = weights(1);
		double w_z = weights(2);

		//std::cout << weights(0) << ", " <<  weights(1) << ", " << weights(2) << std::endl;

		tripletList.push_back(T(i, x + nx*(y + z * ny), (1 - w_x) * (1 - w_y) * (1 - w_z)));
		tripletList.push_back(T(i, x+1 + nx*(y + z * ny), w_x * (1 - w_y) * (1 - w_z)));
		tripletList.push_back(T(i, x + nx*(y+1 + z * ny), (1 - w_x) * w_y * (1 - w_z)));
		tripletList.push_back(T(i, x + nx*(y + (z+1) * ny), (1 - w_x) * (1 - w_y) * w_z));
		tripletList.push_back(T(i, x+1 + nx*(y + (z+1) * ny), w_x * (1 - w_y) * w_z));
		tripletList.push_back(T(i, x+1 + nx*(y+1 + z * ny), w_x * w_y * (1 - w_z)));
		tripletList.push_back(T(i, x + nx*(y+1 + (z+1) * ny), (1 - w_x) * w_y * w_z));
		tripletList.push_back(T(i, x+1 + nx*(y+1 + (z+1) * ny), w_x * w_y * w_z));
	}

	W.setFromTriplets(tripletList.begin(), tripletList.end());
}
