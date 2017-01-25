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
	double xd, yd, zd, xdS, ydS, zdS;
	int Gx, Gy, Gz, GxS, GyS, GzS;
	Eigen::RowVector3d posGrid;
	std::vector<Eigen::Triplet<double>> wVal;
	int sizeWx = (nx - 1)*ny*nz, sizeWy = nx*(ny - 1)*nz, sizeWz = nx*ny*(nz - 1);

	W.resize(numPts, sizeWx + sizeWy + sizeWz);

	for (int i = 0; i < numPts; ++i)
	{
		posGrid = (P.row(i) - corner) / h;
		Gx = floor(posGrid.x());
		Gy = floor(posGrid.y());
		Gz = floor(posGrid.z());
		xd = posGrid.x() - floor(posGrid.x());
		yd = posGrid.y() - floor(posGrid.y());
		zd = posGrid.z() - floor(posGrid.z());

		GxS = convertToStaggered(Gx, xd);
		GyS = convertToStaggered(Gy, yd);
		GzS = convertToStaggered(Gz, zd);
		xdS = convertToStaggered(xd);
		ydS = convertToStaggered(yd);
		zdS = convertToStaggered(zd);

		//construct Wx
		wVal.push_back({ i, (GxS) + (nx - 1)*((Gy) + ny*(Gz)), (1 - xdS)*(1 - yd)*(1 - zd) });
		wVal.push_back({ i, (GxS) + (nx - 1)*((Gy) + ny*(Gz+1)), (1 - xdS)*(1 - yd)*(zd) });
		wVal.push_back({ i, (GxS) + (nx - 1)*((Gy+1) + ny*(Gz)), (1 - xdS)*(yd)*(1 - zd) });
		wVal.push_back({ i, (GxS) + (nx - 1)*((Gy+1) + ny*(Gz+1)), (1 - xdS)*(yd)*(zd) });
		wVal.push_back({ i, (GxS + 1) + (nx - 1)*((Gy) + ny*(Gz)), (xdS)*(1 - yd)*(1 - zd) });
		wVal.push_back({ i, (GxS + 1) + (nx - 1)*((Gy) + ny*(Gz+1)), (xdS)*(1 - yd)*(zd) });
		wVal.push_back({ i, (GxS + 1) + (nx - 1)*((Gy+1) + ny*(Gz)), (xdS)*(yd)*(1 - zd) });
		wVal.push_back({ i, (GxS + 1) + (nx - 1)*((Gy+1) + ny*(Gz+1)), (xdS)*(yd)*(zd) });

		//construct Wy
		wVal.push_back({ i, sizeWx + (Gx) + nx*((GyS) + (ny-1)*(Gz)), (1 - xd)*(1 - ydS)*(1 - zd) });
		wVal.push_back({ i, sizeWx + (Gx) + nx*((GyS) + (ny - 1)*(Gz + 1)), (1 - xd)*(1 - ydS)*(zd) });
		wVal.push_back({ i, sizeWx + (Gx) + nx*((GyS + 1) + (ny - 1)*(Gz)), (1 - xd)*(ydS)*(1 - zd) });
		wVal.push_back({ i, sizeWx + (Gx) + nx*((GyS + 1) + (ny - 1)*(Gz + 1)), (1 - xd)*(ydS)*(zd) });
		wVal.push_back({ i, sizeWx + (Gx + 1) + nx*((GyS) + (ny - 1)*(Gz)), (xd)*(1 - ydS)*(1 - zd) });
		wVal.push_back({ i, sizeWx + (Gx + 1) + nx*((GyS) + (ny - 1)*(Gz + 1)), (xd)*(1 - ydS)*(zd) });
		wVal.push_back({ i, sizeWx + (Gx + 1) + nx*((GyS + 1) + (ny - 1)*(Gz)), (xd)*(ydS)*(1 - zd) });
		wVal.push_back({ i, sizeWx + (Gx + 1) + nx*((GyS + 1) + (ny - 1)*(Gz + 1)), (xd)*(ydS)*(zd) });

		//construct Wz
		wVal.push_back({ i, sizeWx + sizeWy + (Gx) + nx*((Gy) + ny*(GzS)), (1 - xd)*(1 - yd)*(1 - zdS) });
		wVal.push_back({ i, sizeWx + sizeWy + (Gx) + nx*((Gy) + ny*(GzS + 1)), (1 - xd)*(1 - yd)*(zdS) });
		wVal.push_back({ i, sizeWx + sizeWy + (Gx) + nx*((Gy + 1) + ny*(GzS)), (1 - xd)*(yd)*(1 - zdS) });
		wVal.push_back({ i, sizeWx + sizeWy + (Gx) + nx*((Gy + 1) + ny*(GzS + 1)), (1 - xd)*(yd)*(zdS) });
		wVal.push_back({ i, sizeWx + sizeWy + (Gx + 1) + nx*((Gy) + ny*(GzS)), (xd)*(1 - yd)*(1 - zdS) });
		wVal.push_back({ i, sizeWx + sizeWy + (Gx + 1) + nx*((Gy) + ny*(GzS + 1)), (xd)*(1 - yd)*(zdS) });
		wVal.push_back({ i, sizeWx + sizeWy + (Gx + 1) + nx*((Gy + 1) + ny*(GzS)), (xd)*(yd)*(1 - zdS) });
		wVal.push_back({ i, sizeWx + sizeWy + (Gx + 1) + nx*((Gy + 1) + ny*(GzS + 1)), (xd)*(yd)*(zdS) });
	}

	W.setFromTriplets(wVal.begin(), wVal.end());
	W.makeCompressed();
}
