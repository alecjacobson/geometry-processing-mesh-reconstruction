#include "fd_interpolate.h"
#include <iostream>


int index_3d_array1(int i, int j, int k, int nx, int ny)
{
  // An helper function to index W easily. On the assignment note:
  // "We will assume that g_{i,j,k} refers to g(i + j * n_x + k * n_y * n_x)."
  return i + j * nx + k * ny * nx;
}

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  int totalPoint = P.rows();
  // Construct a sparse matrix W with size 'n x (nx * ny * nz)'
  W.resize(totalPoint, nx * ny * nz);

  // Define a tripleList to insert non-zero elements to W.
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(totalPoint * 8);

  for (int i = 0; i < totalPoint; ++i)
  {
      // Get the corresponding grid with corners.
      double x = (P(i, 0) - corner(0)) / h;
      double y = (P(i, 1) - corner(1)) / h;
      double z = (P(i, 2) - corner(2)) / h;
      // Obtain the floor of the corresponding corners.
      int floorX = std::floor(x);
      int floorY = std::floor(y);
      int floorZ = std::floor(z);

      double differenceX = x - floorX;
      double differenceY = y - floorY;
      double differenceZ = z - floorZ;

      // 1. Submit position (floorX, floorY, floorZ).
      tripletList.push_back(T(i, index_3d_array1(floorX, floorY, floorZ, nx, ny),
                             (1 - differenceX) * (1 - differenceY) * (1 - differenceZ)));
      // 2. Submit position (floorX + 1, floorT, floorZ).
      tripletList.push_back(T(i, index_3d_array1(floorX + 1, floorY, floorZ, nx, ny),
                             (differenceX) * (1 - differenceY) * (1 - differenceZ)));
      // 3. Submit position (floorX, floorY + 1, floorZ).
      tripletList.push_back(T(i, index_3d_array1(floorX, floorY + 1, floorZ, nx, ny),
                              (1 - differenceX) * (differenceY) * (1 - differenceZ)));
      // 4. Submit position (floorX, floorY, floorZ + 1).
      tripletList.push_back(T(i, index_3d_array1(floorX, floorY, floorZ + 1, nx, ny),
                              (1 - differenceX) * (1 - differenceY) * (differenceZ)));
      // 5. Submit position (floorX + 1, floorY + 1, floorZ).
      tripletList.push_back(T(i, index_3d_array1(floorX + 1, floorY + 1, floorZ, nx, ny),
                              (differenceX) * (differenceY) * (1 - differenceZ)));
      // 6. Submit position (floorX + 1, floorY, floorZ + 1).
      tripletList.push_back(T(i, index_3d_array1(floorX + 1, floorY, floorZ + 1, nx, ny),
                              (differenceX) * (1 - differenceY) * (differenceZ)));
      // 7. Submit position (floorX, floorY + 1, floorZ + 1).
      tripletList.push_back(T(i, index_3d_array1(floorX, floorY + 1, floorZ + 1, nx, ny),
                              (1 - differenceX) * (differenceY) * (differenceZ)));
      // 8. Submit position (floorX + 1, floorY + 1, floorZ + 1).
      tripletList.push_back(T(i, index_3d_array1(floorX + 1, floorY + 1, floorZ + 1, nx, ny),
                              (differenceX) * (differenceY) * (differenceZ)));
  }
  W.setFromTriplets(tripletList.begin(), tripletList.end());
}
