#include "fd_interpolate.h"
#include <math.h>

using namespace Eigen;

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

	auto getIndxAndWeights = [&](RowVector3d& point)
	{
		RowVector3d indxtemp = (1 / h)*(point - corner);
		RowVector3d indx(floor(indxtemp.x()), floor(indxtemp.y()), floor(indxtemp.z()));
		RowVector3d betaValues = indxtemp - indx;
		RowVectorXd result(6);
		result << indx(0), indx(1), indx(2), betaValues(0), betaValues(1), betaValues(2);
		return result;
	};

	auto makeTriplet = [&](int rowIndx, int i, int j, int k, double value)
	{
		return Triplet<double>(rowIndx, i + j*nx + k*nx*ny, value);
	};

	std::vector<Triplet<double>> tripletList;
	tripletList.reserve(P.rows() * 8);
	W.resize(P.rows(), nx*ny*nz);

	// Reference http://paulbourke.net/miscellaneous/interpolation/
	for (int rowIndx = 0; rowIndx < P.rows(); ++rowIndx)
	{
		RowVector3d point = P.row(rowIndx);
		RowVectorXd info = getIndxAndWeights(point);

		int i = info(0), j = info(1), k = info(2);
		double x = info(3), y = info(4), z = info(5);

		tripletList.push_back(makeTriplet(rowIndx, i, j, k, (1 - x)*(1 - y)*(1 - z)));
		tripletList.push_back(makeTriplet(rowIndx, i + 1, j, k, x*(1 - y)*(1 - z)));
		tripletList.push_back(makeTriplet(rowIndx, i, j + 1, k, (1 - x)*y*(1 - z)));
		tripletList.push_back(makeTriplet(rowIndx, i, j, k + 1, (1 - x)*(1 - y)*z));
		tripletList.push_back(makeTriplet(rowIndx, i + 1, j, k + 1, x*(1 - y)*z));
		tripletList.push_back(makeTriplet(rowIndx, i, j + 1, k + 1, (1 - x)*y*z));
		tripletList.push_back(makeTriplet(rowIndx, i + 1, j + 1, k, x*y*(1 - z)));
		tripletList.push_back(makeTriplet(rowIndx, i + 1, j + 1, k + 1, x*y*z));
	}

	W.setFromTriplets(tripletList.begin(), tripletList.end());
}