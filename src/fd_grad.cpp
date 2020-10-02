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
  ////////////////////////////////////////////////////////////////////////////
	typedef Eigen::SparseMatrix<double> SpMat;
	int Dx_rows = (nx - 1)*ny*nz;
	int Dy_rows = nx*(ny - 1)*nz;
	int Dz_rows = nx*ny*(nz - 1);
	G.resize(Dx_rows + Dy_rows + Dz_rows, nx*ny*nz);
	
	// Compute Dx, Dy, Dz
	SpMat Dx, Dy, Dz;
	fd_partial_derivative(nx, ny, nz, h, 0, Dx);
	fd_partial_derivative(nx, ny, nz, h, 1, Dy);
	fd_partial_derivative(nx, ny, nz, h, 2, Dz);
	
	// Combine Dx, Dy, Dz to G
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(G.rows() * 2);
	for (int col = 0; col < G.outerSize(); ++col) {
		for (SpMat::InnerIterator itDx(Dx, col); itDx; ++itDx)
			triplets.push_back(T(itDx.row(), itDx.col(), itDx.value()));
	}
	for (int col = 0; col < G.outerSize(); ++col) {
		for (SpMat::InnerIterator itDy(Dy, col); itDy; ++itDy)
			triplets.push_back(T(Dx_rows + itDy.row(), itDy.col(), itDy.value()));
	}
	for (int col = 0; col < G.outerSize(); ++col){
		for (SpMat::InnerIterator itDz(Dz, col); itDz; ++itDz)
			triplets.push_back(T(Dx_rows + Dy_rows + itDz.row(), itDz.col(), itDz.value()));
	}
	G.setFromTriplets(triplets.begin(), triplets.end());
}
