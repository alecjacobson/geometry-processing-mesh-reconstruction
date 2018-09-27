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

  auto ind = [&](int x, int y, int z) {
    return x + nx * (y + z * ny);
  };

  // Calculate m
  int m;
  switch (dir) {
    case 0:
      m = (nx - 1) * ny * nz;
      break;
    case 1:
      m = nx * (ny - 1) * nz;
      break;
    case 2:
      m = nx * ny * (nz - 1);
      break;
    default:
      m = nx * ny * (nz - 1);
  }

  D.resize(m, nx * ny * nz);

  // Calculate D
  for (int i = 0; i < m; i++) {
    D.insert(i, i) = 1 / h;
    D.insert(i, i + 1) = -1 / h;
  }

  ////////////////////////////////////////////////////////////////////////////
}
