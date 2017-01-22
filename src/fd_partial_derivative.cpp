#include "fd_partial_derivative.h"

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
	tripletList.reserve(nx*ny*nz);

	int iStart = 0;
	int jStart = 0;
	int kStart = 0;
	switch (dir) {
		case 0: iStart = 1; break;
		case 1: jStart = 1; break;
		case 3: kStart = 1; break;
		default: break;
	}

	for (int i = iStart; i < nx; ++ i) {
		for (int j = jStart; j < ny; ++j) {
			for (int k = kStart; k < nz; ++k) {
				int row = i + j*nx + k*nx*ny;
				for (int l = 0; l < nx*ny*nz; ++l) {
					int dirIndex = 0;
					switch (dir) {
						case 0:  dirIndex = i; break;
						case 1:  dirIndex = j; break;
						case 2:  dirIndex = k; break;
						default: dirIndex = i; break;
					}
					
					if (l == dirIndex - 1) {
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
