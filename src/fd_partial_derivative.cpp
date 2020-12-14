#include "fd_partial_derivative.h"

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
	
	int m;
	if (dir == 0) {
		m = (nx - 1)*ny*nz;
	}
	else if (dir == 1) {
		m = nx * (ny - 1)*nz;
	}
	else {
		m = nx * ny*(nz - 1);
	}
	Eigen::SparseMatrix<double> Dx(m, nx*ny*nz);
	Dx.setZero();
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(m*2);

	// Form the finite difference matrices based
	int curr_row = 0;
	int curr_col = 0;
	if (dir == 0) {
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					if (i == nx - 1) {
						curr_col += 1;
					}
					else {
						triplets.push_back({ curr_row, curr_col  , -1  });
						triplets.push_back({ curr_row, curr_col + 1,  1  });
						curr_row++; curr_col++;
					}
				}
			}
		}
	}
	else if (dir == 1) {
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					if (j == ny - 1) {
						curr_col += 1;
					}	else {
						triplets.push_back({ curr_row, curr_col  , -1 });
						triplets.push_back({ curr_row, curr_col + nx,  1  });
						curr_row++; curr_col++;
					}
				}
			}
		}
	}
	else {
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					if (k == nz - 1) {
						curr_col += 1;
					}
					else {

						triplets.push_back({ curr_row, i + j * nx + k * nx*ny  , -1  });
						triplets.push_back({ curr_row, i + j * nx + k * nx*ny + nx * ny,  1  });
						curr_row++;
					}
				}
			}
		}
	}
	
	Dx.setFromTriplets(triplets.begin(), triplets.end());
	D = Dx;
  ////////////////////////////////////////////////////////////////////////////
}
