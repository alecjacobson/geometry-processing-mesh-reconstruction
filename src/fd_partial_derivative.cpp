#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
	int NX = nx, NY = ny, NZ = nz;
	std::vector<Eigen::Triplet<double>> dVal;

	if (dir == 1)
		NX--;
	else if (dir == 2)
		NY--;
	else
		NZ--;

	D.resize(NX*NY*NZ, nx*ny*nz);

	
}
