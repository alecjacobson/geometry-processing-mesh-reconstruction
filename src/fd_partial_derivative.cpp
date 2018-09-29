#include "fd_partial_derivative.h"

int index_3d_array(int i, int j, int k, int nx, int ny)
{
    // We will assume that g_{i,j,k} refers to g(i + j * n_x + k * n_y * n_x)
    return i + j * nx + k * ny * nx;
}

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  // m: number of rows in D.
  int m;
  int modifiedNx = nx;
  int modifiedNy = ny;
  int modifiedNz = nz;
  // dir should be 0, 1 or 2.
  if (dir == 0)
  {
      m = (nx - 1) * ny * nz;
      modifiedNx = (nx - 1);
  } else if (dir == 1)
  {
      m = nx * (ny - 1) * nz;
      modifiedNy = (ny - 1);
  } else
  {
      m = nx * ny * (nz - 1);
      modifiedNz = (nz - 1);
  }
  // From assignment note:
  //  D: m by nx*ny*nz sparse partial derivative matrix, where:
  //     m = (nx-1)*ny*nz  if dir = 0
  //     m = nx*(ny-1)*nz  if dir = 1
  //     m = nx*ny*(nz-1)  otherwise (if dir = 2)
  D.resize(m, nx * ny * nz);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  int row, col;
  int leftIndex, rightIndex;
  for (int i = 0; i < modifiedNx; ++i)
  {
      for (int j = 0; j < modifiedNy; ++j)
      {
          for (int k = 0; k < modifiedNz; ++k)
          {
              row = index_3d_array(i, j, k, modifiedNx, modifiedNy);
              tripletList.push_back(T(row, index_3d_array(i, j, k, nx, ny), -1.0));
              if (dir == 0)
              {
                  col = index_3d_array(i + 1, j, k, nx, ny);
              } else if (dir == 1)
              {
                  col = index_3d_array(i, j + 1, k, nx, ny);
              } else
              {
                  col = index_3d_array(i, j, k + 1, nx, ny);
              }
              tripletList.push_back(T(row, col, 1.0));
          }
      }
  }
  D.setFromTriplets(tripletList.begin(), tripletList.end());
}
