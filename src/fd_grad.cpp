#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <igl/cat.h>

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

    // matrix: G = [Dx;Dy;Dz]
	int m = (nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1);
	G.resize(m, nx*ny*nz);

	Eigen::SparseMatrix<double> Dx;
	Eigen::SparseMatrix<double> Dy;
	Eigen::SparseMatrix<double> Dz;

	fd_partial_derivative(nx, ny, nz, h, 0, Dx);
	fd_partial_derivative(nx, ny, nz, h, 1, Dy);
	fd_partial_derivative(nx, ny, nz, h, 2, Dz);

	G = igl::cat(1, Dx, Dy);
	G = igl::cat(1, G, Dz);

}
