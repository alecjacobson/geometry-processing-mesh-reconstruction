#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
	const int nx,
	const int ny,
	const int nz,
	const double h,
	Eigen::SparseMatrix<double> & G)
{
	////////////////////////////////////////////////////////////////////////////
	// Add your code here
	Eigen::SparseMatrix<double> Gx, Gy, Gz;

	// Calculate partial derivative matrix for each x, y, z direction
	fd_partial_derivative(nx, ny, nz, h, 0, Gx);
	fd_partial_derivative(nx, ny, nz, h, 1, Gy);
	fd_partial_derivative(nx, ny, nz, h, 2, Gz);

	// Concatenate the matrix vertically
	Eigen::SparseMatrix<double> M(Gx.rows() + Gy.rows() + Gz.rows(), Gx.cols());
	M.reserve(Gx.nonZeros() + Gy.nonZeros() + Gz.nonZeros());
	for (int c = 0; c < Gx.cols(); c++) {
		M.startVec(c);
		for (Eigen::SparseMatrix<double>::InnerIterator itx(Gx, c); itx; ++itx) {
			M.insertBack(itx.row(), c) = itx.value();
		}
		for (Eigen::SparseMatrix<double>::InnerIterator ity(Gy, c); ity; ++ity) {
			M.insertBack(ity.row() + Gx.rows(), c) = ity.value();
		}
		for (Eigen::SparseMatrix<double>::InnerIterator itz(Gz, c); itz; ++itz) {
			M.insertBack(itz.row() + Gx.rows() + Gy.rows(), c) = itz.value();
		}		
	}
	////////////////////////////////////////////////////////////////////////////
	G = M;
	M.finalize();
}
