#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
	int numRow;
	int numCol = nx*ny*nz;
	int sizeDx, sizeDy;
	Eigen::SparseMatrix<double> Dx, Dy, Dz;
	std::vector<Eigen::Triplet<double>> gVal;

	fd_partial_derivative(nx, ny, nz, h, 1, Dx);
	fd_partial_derivative(nx, ny, nz, h, 2, Dy);
	fd_partial_derivative(nx, ny, nz, h, 3, Dz);

	sizeDx = Dx.rows();
	sizeDy = Dy.rows();
	numRow = sizeDx + sizeDy + Dz.rows();

	G.resize(numRow, numCol);
	G.makeCompressed();
	G.reserve(2 * numRow);

	for (int k = 0; k < Dx.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(Dx, k); it; ++it)
			gVal.push_back({it.row(), it.col(), it.value()});

	for (int k = 0; k < Dy.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(Dy, k); it; ++it)
			gVal.push_back({ sizeDx + it.row(), it.col(), it.value() });

	for (int k = 0; k < Dz.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(Dz, k); it; ++it)
			gVal.push_back({ sizeDx + sizeDy + it.row(), it.col(), it.value() });
	
	G.setFromTriplets(gVal.begin(), gVal.end());
}
