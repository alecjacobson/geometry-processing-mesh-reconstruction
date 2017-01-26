#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
	assert(dir >= 1 && dir <= 3);

	int NX = nx, NY = ny, NZ = nz;
	std::vector<Eigen::Triplet<double>> dVal;

	if (dir == 1)
		NX--;
	else if (dir == 2)
		NY--;
	else if (dir == 3)
		NZ--;

	D.resize(NX*NY*NZ, nx*ny*nz);

	for (int i = 0; i < NX; ++i)
		for(int j=0; j < NY; ++j)
			for (int k = 0; k < NZ; ++k)
			{
				switch (dir)
				{
					//x
				case 1:
					dVal.push_back({i + NX*(j + k*NY), i + 1 + nx*(j + k*ny), 1 });
					dVal.push_back({i + NX*(j + k*NY), i + nx*(j + k*ny), -1});
					break;
					//y
				case 2:
					dVal.push_back({ i + NX*(j + k*NY), i + nx*(j + 1 + k*ny), 1 });
					dVal.push_back({ i + NX*(j + k*NY), i + nx*(j + k*ny), -1 });
					break;
					//z
				case 3:
					dVal.push_back({ i + NX*(j + k*NY), i + nx*(j + (k + 1)*ny), 1 });
					dVal.push_back({ i + NX*(j + k*NY), i + nx*(j + k*ny), -1 });
					break;
				}
			}

	D.setFromTriplets(dVal.begin(), dVal.end());
}
