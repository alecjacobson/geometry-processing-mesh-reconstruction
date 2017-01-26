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
	tripletList.reserve(2*D.rows());

	int iEnd = nx;
	int jEnd = ny;
	int kEnd = nz;

	switch (dir) {
		case 0: iEnd--; break;
		case 1: jEnd--; break;
		case 2: kEnd--; break;
	}

	for (int i = 0; i < iEnd; ++i) {
		for (int j = 0; j < jEnd; ++j) {
			for (int k = 0; k < kEnd; ++k) {
				// Not too sure about these column indicies
				switch (dir) {
					case 0: 
						tripletList.push_back(Eigen::Triplet<double>(i + j*iEnd + k*iEnd*ny, i + j*nx + k*nx*ny, -1/h));
						tripletList.push_back(Eigen::Triplet<double>(i + j*iEnd + k*iEnd*ny, i + 1 + j*nx + k*nx*ny, 1/h));
						break;
					case 1:
						tripletList.push_back(Eigen::Triplet<double>(i + j*nx + k*nx*jEnd, i + j*nx + k*nx*ny, -1/h));
						tripletList.push_back(Eigen::Triplet<double>(i + j*nx + k*nx*jEnd, i + (j+1)*nx + k*nx*ny, 1/h));
						break;
					case 2:
						tripletList.push_back(Eigen::Triplet<double>(i + j*nx + k*nx*ny, i + j*nx + k*nx*ny, -1/h));
						tripletList.push_back(Eigen::Triplet<double>(i + j*nx + k*nx*ny, i + j*nx + (k + 1)*nx*ny, 1/h));
						break;
				}
				
			}
		}
	}

	D.setFromTriplets(tripletList.begin(), tripletList.end());
	
}
