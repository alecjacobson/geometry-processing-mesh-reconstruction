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
	const int n = P.rows();

	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(n);
	
	for (int t = 0; t < n; ++t) {
		double x = P(t, 0);
		double y = P(t, 1);
		double z = P(t, 2);

		int i = floor((x - corner(0)) / h);
		int j = floor((y - corner(1)) / h);
		int k = floor((z - corner(2)) / h);

		double x0 = corner(0) + i*h;
		double x1 = corner(0) + (i+1)*h;

		double y0 = corner(1) + j*h;
		double y1 = corner(1) + (j+1)*h;

		double z0 = corner(2) + k*h;
		double z1 = corner(2) + (k+1)*h;

		double xd = (x - x0) / (x1 - x0);
		double yd = (y - y0) / (y1 - y0);
		double zd = (z - z0) / (z1 - z0);

		/*
		Eigen::RowVector3d Vx0y0z0 = grid.row(i + j*nx + k*nx*ny);
		Eigen::RowVector3d Vx1y0z0 = grid.row(i + 1 + j*nx + k*nx*ny);
		Eigen::RowVector3d Vx0y0z1 = grid.row(i + j*nx + (k + 1)*nx*ny);
		Eigen::RowVector3d Vx1y0z1 = grid.row(i + 1 + j*nx + (k + 1)*nx*ny);
		Eigen::RowVector3d Vx0y1z0 = grid.row(i + (j + 1)*nx + k*nx*ny);
		Eigen::RowVector3d Vx1y1z0 = grid.row(i + 1 + (j + 1)*nx + k*nx*ny);
		Eigen::RowVector3d Vx0y1z1 = grid.row(i + (j + 1)*nx + (k + 1)*nx*ny);
		Eigen::RowVector3d Vx1y1z1 = grid.row(i + 1 + (j + 1)*nx + (k + 1)*nx*ny);

		Eigen::RowVector3d c00 = Vx0y0z0*(1 - xd) + Vx1y0z0*xd;
		Eigen::RowVector3d c01 = Vx0y0z1*(1 - xd) + Vx1y0z1*xd;
		Eigen::RowVector3d c10 = Vx0y1z0*(1 - xd) + Vx1y1z0*xd;
		Eigen::RowVector3d c11 = Vx0y1z1*(1 - xd) + Vx1y1z1*xd;

		Eigen::RowVector3d c0 = c00*(1 - yd) + c10*yd;
		Eigen::RowVector3d c1 = c01*(1 - yd) + c11*yd;

		Eigen::RowVector3d c = c0*(1 - zd) + c1*zd;

		Pattern: x0 -> (1-xd) in weight, x1 -> xd in weight. Same for y and z.

		*/
		
		tripletList.push_back(Eigen::Triplet<double>(t, i + j*nx + k*nx*ny, (1 - xd)*(1 - yd)*(1 - zd))); //(x0,y0,z0)
		tripletList.push_back(Eigen::Triplet<double>(t, i + 1 + j*nx + k*nx*ny, (xd)*(1 - yd)*(1 - zd))); //(x1,y0,z0)

		tripletList.push_back(Eigen::Triplet<double>(t, i + j*nx + (k + 1)*nx*ny, (1 - xd)*(1 - yd)*(zd))); //(x0,y0,z1)
		tripletList.push_back(Eigen::Triplet<double>(t, i + 1 + j*nx + k*x*ny, (xd)*(1 - yd)*(zd))); //(x1,y0,z1)

		tripletList.push_back(Eigen::Triplet<double>(t, i + (j + 1)*nx + k*nx*ny, (1 - xd)*(yd)*(1 - zd))); //(x0,y1,z0)
		tripletList.push_back(Eigen::Triplet<double>(t, i + 1 + (j + 1)*nx + k*x*ny, (xd)*(yd)*(1 - zd))); //(x1,y1,z0)

		tripletList.push_back(Eigen::Triplet<double>(t, i + (j + 1)*nx + (k + 1)*nx*ny, (1 - xd)*(yd)*(zd))); //(x0,y1,z1)
		tripletList.push_back(Eigen::Triplet<double>(t, i + 1 + (j + 1)*nx + (k + 1)*nx*ny, (xd)*(yd)*(zd))); //(x1,y1,z1)
		
	}

	W.setFromTriplets(tripletList.begin(), tripletList.end());
	W.makeCompressed();
}
