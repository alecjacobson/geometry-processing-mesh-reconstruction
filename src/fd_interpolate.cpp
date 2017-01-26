#include "fd_interpolate.h"
#include <cmath>

inline double convertToStaggered(double x)
{
	return (x < 0.5) ? x + 0.5 : x - 0.5;
}

inline int convertToStaggered(int Gx, double frac)
{
	return (frac < 0.5) ? Gx : Gx + 1;
}

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
	int numPts = P.rows();
	double xd, yd, zd;
	int Gx, Gy, Gz;
	Eigen::RowVector3d posGrid;
	std::vector<Eigen::Triplet<double>> wVal;

	W.resize(numPts, nx*ny*nz);
	
	for (int i = 0; i < numPts; ++i)
	{
		posGrid = (P.row(i) - corner) / h;

		Gx = floor(posGrid.x());
		Gy = floor(posGrid.y());
		Gz = floor(posGrid.z());

		xd = posGrid.x() - Gx;
		yd = posGrid.y() - Gy;
		zd = posGrid.z() - Gz;

		wVal.push_back({ i, Gx + nx*(Gy + ny*Gz), (1 - xd)*(1 - yd)*(1 - zd) });
		wVal.push_back({ i, Gx + nx*(Gy + ny*(Gz+1)), (1 - xd)*(1 - yd)*(zd) });
		wVal.push_back({ i, Gx + nx*(Gy+1 + ny*Gz), (1 - xd)*(yd)*(1 - zd) });
		wVal.push_back({ i, Gx + nx*(Gy+1 + ny*(Gz+1)), (1 - xd)*(yd)*(zd) });
		wVal.push_back({ i, Gx+1 + nx*(Gy + ny*Gz), (xd)*(1 - yd)*(1 - zd) });
		wVal.push_back({ i, Gx+1 + nx*(Gy + ny*(Gz+1)), (xd)*(1 - yd)*(zd) });
		wVal.push_back({ i, Gx+1 + nx*(Gy+1 + ny*Gz), (xd)*(yd)*(1 - zd) });
		wVal.push_back({ i, Gx+1 + nx*(Gy+1 + ny*(Gz+1)), (xd)*(yd)*(zd) });
	}

	W.setFromTriplets(wVal.begin(), wVal.end());
}
