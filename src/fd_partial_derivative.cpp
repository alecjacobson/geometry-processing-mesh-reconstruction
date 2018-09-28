#include "fd_partial_derivative.h"
typedef Eigen::Triplet<double> T;

void fd_partial_derivative(
    const int nx,
    const int ny,
    const int nz,
    const double h,
    const int dir,
    Eigen::SparseMatrix<double> &D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  int m = 0;
  int space = 0;

  switch (dir)
  {
  case 0:
    m = (nx - 1) * ny * nz;
    space = 1;
    break;

  case 1:
    m = nx * (ny - 1) * nz;
    space = nx;
    break;

  case 2:
    m = nx * ny * (nz - 1);
    space = nx * ny;
    break;
  }

  D.resize(m, nx * ny * nz);

  std::vector<T> tlist;
  tlist.reserve(2 * m);

  int count = 0;

  for (int k = 0; k < nz; k++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {
        int index = i + j * nx + k * ny * nx;
        int col = index;
        int row = count;

        if ((i == nx - 1 && dir == 0) || (j == ny - 1 && dir == 1) || (k == nz - 1 && dir == 2))
        {
          //skip where index reach last point in a loop
          continue;
        }
        else
        {
          tlist.push_back(T(row, col, -1.0 / h));
          tlist.push_back(T(row, col + space, 1.0 / h));
          count++;
        }
      }
    }
  }

  D.setFromTriplets(tlist.begin(), tlist.end());
}
