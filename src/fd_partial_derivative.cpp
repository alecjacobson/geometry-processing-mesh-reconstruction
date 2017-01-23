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
	//O((nx*ny*nz)^2)
	std::vector<Eigen::Triplet<double>> tripletList;
	tripletList.reserve(D.rows());

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			for (int k = 0; k < nz; ++k) {
				int row = row = i + j*nx+ k*nx*ny;
				int dirIndex = 0;
				switch (dir) {
					case 0:  dirIndex = i; break;
					case 1:  dirIndex = j; break;
					case 2:  dirIndex = k; break;
					default: dirIndex = i; break;
				}

				if (dirIndex > 0) {
					tripletList.push_back(Eigen::Triplet<double>(row, dirIndex - 1, -1 / h));
				}
				tripletList.push_back(Eigen::Triplet<double>(row, dirIndex, 1/h));
				
			}
		}
	}

	D.setFromTriplets(tripletList.begin(), tripletList.end());
	D.makeCompressed();
	
}
