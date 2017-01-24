#include "fd_grad.h"

#include "fd_partial_derivative.h"
void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
	//We simply use fd_partial_derivatives to compute the submatrices, 
	//and concatenate them on the fly.
	const int rows = (nx - 1)*ny*nz + nx*(ny - 1)*nz + nx*ny*(nz - 1);
	const int cols = nx*ny*nz;
	Eigen::SparseMatrix<double> D;
	std::vector<Eigen::Triplet<double>> entries;
	entries.reserve(2*rows);
	int row_offset = 0;
	for (int i = 0; i < 3; i++) {
		Eigen::RowVector3d offset = Eigen::RowVector3d::Zero();
		offset(i) = 1;

		fd_partial_derivative(nx, ny, nz, h, i, D);
		for (int k = 0; k < D.outerSize(); k++) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(D, k); it; ++it) {
				entries.emplace_back(it.row() + row_offset, it.col(), it.value());
			}
		}
		row_offset += D.rows();
	}
	
	G.resize(rows, cols);
	G.setFromTriplets(entries.cbegin(), entries.cend());
}
