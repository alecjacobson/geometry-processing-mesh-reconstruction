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

  int nxm = nx;
  int nym = ny;
  int nzm = nz;

  int m;
  switch (dir) {
    case 0:
      nxm = nx - 1;
      break;
    case 1:
      nym = ny - 1;
      break;
    default:
      nzm = nz - 1;
  }
  m = nxm * nym * nzm;

  auto ind = [&](int x, int y, int z) {
      return x + nx * (y + z * ny);
  };
  auto mind = [&](int x, int y, int z) {
      return x + nxm * (y + z * nym);
  };

  D.resize(m, nx * ny * nz);

  // Calculate D
  for (int i = 0; i < nxm; i++) {
    for (int j = 0; j < nym; j++) {
      for (int k = 0; k < nzm; k++) {
        D.insert(mind(i, j, k), ind(i, j, k)) = - 1 / h;
        switch(dir) {
          case 0:
            D.insert(mind(i, j, k), ind(i + 1, j, k)) = 1 / h;
            break;
          case 1:
            D.insert(mind(i, j, k), ind(i, j + 1, k)) = 1 / h;
            break;
          default:
            D.insert(mind(i, j, k), ind(i, j, k + 1)) = 1 / h;
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////
}
