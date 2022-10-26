#include "fd_partial_derivative.h"
#include <iostream>

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

	int m = 0;
	int n_nx = nx;
	int n_ny = ny;
	int n_nz = nz;

	switch (dir)
	{
	case 0:
		m = (nx - 1)*ny*nz;//  if dir = 0
		n_nx -= 1;
		break;
	case 1:
		m = nx*(ny - 1)*nz;//  if dir = 1
		n_ny -= 1;
		break;
	case 2:
		m = nx*ny*(nz - 1);//  if dir = 2
		n_nz -= 1;
		break;
	default:
		m = 0;
		break;
	}

	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(m);

	D.resize(m, nx*ny*nz);

	for (int i = 0; i < n_nx; i++)
	{
		for (int j = 0; j < n_ny; j++)
		{
			for (int k = 0; k < n_nz; k++)
			{
				int l = i + j*n_nx + k*n_ny*n_nx;
				int t = i + j*nx + k*ny*nx;

				if (dir == 0)
				{
					tripletList.push_back(T(l, t, -1./h));
					tripletList.push_back(T(l, (i + 1) + j*nx + k*ny*nx, 1./h));
				}
				else if (dir == 1)
				{
					tripletList.push_back(T(l, t, -1. / h));
					tripletList.push_back(T(l, i + (j + 1)*nx + k*ny*nx, 1. / h));
				}
				else
				{
					tripletList.push_back(T(l, t, -1. / h));
					tripletList.push_back(T(l, i + j*nx + (k + 1)*ny*nx, 1. / h));
				}
			}
		}
	}

	D.setFromTriplets(tripletList.begin(), tripletList.end());

}
