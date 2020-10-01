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
  
  // Compute the derivatives for the desired direction
  switch (dir)
  {
  	case 1:
  	{
  		D.resize((nx-1)*ny*nz, nx*ny*nz);

		// Using the triplet insertion method based on Eigen's documentation
		typedef Eigen::Triplet<double> T;
		std::vector<T> tripletList;

		for (int ii = 0; ii < nx*ny*nz; ii++)
		{
			int col = ii;
			int col_other;

			// Check if we're at an edge, in which case the finite difference is taken to be a forward difference
			if (ii%nx == 0)
			{
				col_other = ii+1;
				if (col >= (nx-1)*ny*nz)
					col = (nx-1)*ny*nz - 1;
				if (col_other >= nx*ny*nz)
					col_other = nx*ny*nz - 1;

				tripletList.push_back(T(col, col, -1.0));
				tripletList.push_back(T(col, col_other, 1.0));
			}
			else
			{
				col_other = ii-1;
				if (col_other >= (nx-1)*ny*nz)
					col_other = (nx-1)*ny*nz - 1;
				if (col_other < 0)
					col_other = 0;

				tripletList.push_back(T(col_other, col_other, -1.0));
				tripletList.push_back(T(col_other, col, 1.0));
			}
		}

  		D.setFromTriplets(tripletList.begin(), tripletList.end());
  		D.finalize();
  	}
  	case 2:
  	{
  		D.resize(nx*(ny-1)*nz, nx*ny*nz);

		// Using the triplet insertion method based on Eigen's documentation
		typedef Eigen::Triplet<double> T;
		std::vector<T> tripletList;

		for (int ii = 0; ii < nx*ny*nz; ii++)
		{
			int col = ii;
			int col_other;

			// Check if we're at an edge, in which case the finite difference is taken to be a forward difference
			if (ii%(nx*ny) < nx)
			{
				col_other = ii+1*nx;
				if (col >= nx*(ny-1)*nz)
					col = nx*(ny-1)*nz - 1;
				if (col_other >= nx*ny*nz)
					col_other = nx*ny*nz - 1;

				tripletList.push_back(T(col, col, -1.0));
				tripletList.push_back(T(col, col_other, 1.0));
			}
			else
			{
				col_other = ii-1*nx;
				if (col_other >= nx*(ny-1)*nz)
					col_other = nx*(ny-1)*nz - 1;
				if (col_other < 0)
					col_other = 0;

				tripletList.push_back(T(col_other, col_other, -1.0));
				tripletList.push_back(T(col_other, col, 1.0));
			}
		}

  		D.setFromTriplets(tripletList.begin(), tripletList.end());
  		D.finalize();
  	}
  	case 3:
  	{
  		D.resize(nx*ny*(nz-1), nx*ny*nz);

		// Using the triplet insertion method based on Eigen's documentation
		typedef Eigen::Triplet<double> T;
		std::vector<T> tripletList;

		for (int ii = 0; ii < nx*ny*nz; ii++)
		{
			int col = ii;
			int col_other;

			// Check if we're at an edge, in which case the finite difference is taken to be a forward difference
			if (ii < nx)
			{
				col_other = ii+1*nx*ny;
				if (col >= nx*ny*(nz-1))
					col = nx*ny*(nz-1) - 1;
				if (col_other >= nx*ny*nz)
					col_other = nx*ny*nz - 1;

				tripletList.push_back(T(col, col, -1.0));
				tripletList.push_back(T(col, col_other, 1.0));
			}
			else
			{
				col_other = ii-1*nx*ny;
				if (col_other >= nx*ny*(nz-1))
					col_other = nx*ny*(nz-1) - 1;
				if (col_other < 0)
					col_other = 0;

				tripletList.push_back(T(col_other, col_other, -1.0));
				tripletList.push_back(T(col_other, col, 1.0));

			}
		}

  		D.setFromTriplets(tripletList.begin(), tripletList.end());  		
  		D.finalize();
  	}

  };


  ////////////////////////////////////////////////////////////////////////////
}
