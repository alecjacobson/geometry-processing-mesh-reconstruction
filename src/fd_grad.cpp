#include "fd_grad.h"
#include "fd_partial_derivative.h"

typedef Eigen::Triplet<double> Triplet;

/**
* @brief Concatenates vertically.
*/
void concatSparse(Eigen::SparseMatrix<double>& source,
		std::vector<Triplet>& targetTriplets, int firstTargetRow) {
	typedef Eigen::SparseMatrix<double>::InnerIterator Iterator;

	// From Eigen tutorials...
	for (int k = 0; k < source.outerSize(); k++) {
		for (Iterator it(source, k); it; ++it) {
			int sourceRow = it.row();
			int column = it.col();
			double value = it.value();
			targetTriplets.push_back(Triplet(firstTargetRow + sourceRow, column, value));
		}
	}

}

void fd_grad(
		const int nx,
		const int ny,
		const int nz,
		const double h,
		Eigen::SparseMatrix<double> & G) {

	Eigen::SparseMatrix<double> Dx((nx - 1) * ny * nz, nx * ny * nz);
	fd_partial_derivative(nx, ny, nz, h, 0, Dx);

	Eigen::SparseMatrix<double> Dy(nx * (ny - 1) * nz, nx * ny * nz);
	fd_partial_derivative(nx, ny, nz, h, 1, Dy);

	Eigen::SparseMatrix<double> Dz(nx * ny * (nz - 1), nx * ny * nz);
	fd_partial_derivative(nx, ny, nz, h, 2, Dz);

	// Concatenate
	std::vector<Triplet> triplets;
	concatSparse(Dx, triplets, 0);
	concatSparse(Dy, triplets, Dx.rows());
	concatSparse(Dz, triplets, Dx.rows() + Dy.rows());	

	G.resize(Dx.rows() + Dy.rows() + Dz.rows(), Dx.cols());
	G.setFromTriplets(triplets.begin(), triplets.end());

}
