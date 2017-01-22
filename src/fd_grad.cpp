#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
	Eigen::SparseMatrix<double> Dx((nx-1)*ny*nz,nx*ny*nz);
	Eigen::SparseMatrix<double> Dy(nx*(ny-1)*nz, nx*ny*nz);
	Eigen::SparseMatrix<double> Dz(nx*ny*(nz-1), nx*ny*nz);

	fd_partial_derivative(nx, ny, nz, h, 0, Dx);
	fd_partial_derivative(nx, ny, nz, h, 1, Dy);
	fd_partial_derivative(nx, ny, nz, h, 2, Dz);

	/* This seems like it should work but it will not compile 
	G.topRows((nx - 1)*ny*nz) = Dx;
	G.middleRows((nx-1)*ny*nz, nx*(ny-1)*nz) = Dy;
	G.bottomRows(nx*ny*(nz-1)) = Dz;
	*/

	for (int c = 0; c < nx*ny*nz; ++c) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(Dx, c); it; ++it) {
			G.insertBack(it.row(), c) = it.value();
		}
		for (Eigen::SparseMatrix<double>::InnerIterator it(Dy, c); it; ++it) {
			G.insertBack(it.row(), c) = it.value();
		}
		for (Eigen::SparseMatrix<double>::InnerIterator it(Dz, c); it; ++it) {
			G.insertBack(it.row(), c) = it.value();
		}
	}
	G.finalize();

}
