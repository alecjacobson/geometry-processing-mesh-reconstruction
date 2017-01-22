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
	tripletList.reserve(2*nx*ny*nz);

	int iEnd = nx;
	int jEnd = ny;
	int kEnd = 0;
	switch (dir) {
		case 0: iEnd--; break;
		case 1: jEnd--; break;
		case 2: kEnd--; break;
		default: break;
	}

	for (int i = 0; i < iEnd; ++ i) {
		for (int j = 0; j < jEnd; ++j) {
			for (int k = 0; k < kEnd; ++k) {
				int row = i + j*nx + k*nx*ny;
				for (int l = 0; l < nx*ny*nz; ++l) {
					int dirIndex = 0;
					switch (dir) {
						case 0:  dirIndex = i; break;
						case 1:  dirIndex = j; break;
						case 2:  dirIndex = k; break;
						default: dirIndex = i; break;
					}

					if (row >= D.rows()) {
						std::cout << row << std::endl;
					}
					
					if (l == (dirIndex - 1)) {
						tripletList.push_back(Eigen::Triplet<double>(row, l, -1/h));
					}
					else if(l == dirIndex) {
						tripletList.push_back(Eigen::Triplet<double>(row, l, 1/h));
					}
				}
			}
		}
	}

	D.setFromTriplets(tripletList.begin(), tripletList.end());
	D.makeCompressed();
	
}
