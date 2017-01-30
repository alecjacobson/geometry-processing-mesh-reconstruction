#include "fd_partial_derivative.h"

using namespace Eigen;

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

	int stag_nx = nx, stag_ny = ny, stag_nz = nz;

	if (dir == 0) stag_nx = nx - 1;
	else if (dir == 1) stag_ny = ny - 1;
	else stag_nz = nz - 1;

	int D_rows = stag_nx*stag_ny*stag_nz;

	std::vector<Triplet<double>> tripletList;
	tripletList.reserve(D_rows * 2);

	auto makeTriplet = [&](int i, int j, int k, int l, double value)
	{
		return Triplet<double>(i + j*stag_nx + k*stag_nx*stag_ny, l, value);
	};

	for (int i = 0; i < stag_nx; ++i)
	{
		for (int j = 0; j < stag_ny; ++j)
		{
			for (int k = 0; k < stag_nz; ++k)
			{
				int l1 = i + j*nx + k*nx*ny;
				int l2;

				if (dir == 0) l2 = (i + 1) + j*nx + k*nx*ny;
				else if (dir == 1) l2 = i + (j + 1)*nx + k*nx*ny;
				else l2 = i + j*nx + (k + 1)*nx*ny;

				tripletList.push_back(makeTriplet(i, j, k, l1, -1 / h));
				tripletList.push_back(makeTriplet(i, j, k, l2, 1 / h));
			}
		}
	}

	D.resize(D_rows, nx*ny*nz);
	D.setFromTriplets(tripletList.begin(), tripletList.end());
}
