#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
	Eigen::SparseMatrix<double> Dx((nx-1)*ny*nz, nx*ny*nz);
	Eigen::SparseMatrix<double> Dy(nx*(ny-1)*nz, nx*ny*nz);
	Eigen::SparseMatrix<double> Dz(nx*ny*(nz-1), nx*ny*nz);

	fd_partial_derivative(nx, ny, nz, h, 0, Dx);
	fd_partial_derivative(nx, ny, nz, h, 1, Dy);
	fd_partial_derivative(nx, ny, nz, h, 2, Dz);

	std::vector<Eigen::Triplet<double> > tripletList;
	//tripletList.reserve(G.rows());
	for (int c = 0; c < Dx.outerSize(); ++c) {
		 for (Eigen::SparseMatrix<double>::InnerIterator it(Dx, c); it; ++it) {
			tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
		 }
	}
	for (int c = 0; c < Dy.outerSize(); ++c) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(Dy, c); it; ++it) {
			tripletList.push_back(Eigen::Triplet<double>((nx - 1)*ny*nz + it.row(), it.col(), it.value()));
		}
	}
	for (int c = 0; c < Dz.outerSize(); ++c) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(Dz, c); it; ++it) {
			tripletList.push_back(Eigen::Triplet<double>((nx - 1)*ny*nz + nx*(ny - 1)*nz + it.row(), it.col(), it.value()));
		}
	}
	G.setFromTriplets(tripletList.begin(), tripletList.end());

}
