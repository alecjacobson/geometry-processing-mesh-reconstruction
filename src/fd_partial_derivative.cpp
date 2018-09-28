#include "fd_partial_derivative.h"
#include <vector>

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir, //0 = x, 1 = y, 2 = z
  Eigen::SparseMatrix<double> & D)
{

  int x = nx;
  int y = ny;
  int z = nz;
  if (dir == 0) {
    x--;
  } else if (dir == 1) {
    y--;
  } else if (dir == 2) {
    z--;
  }

  int ind;
  std::vector<Eigen::Triplet<double>> triplets; 
  for(int i = 0; i < x; i++) 
  {
    for(int j = 0; j < y; j++)
    {
      for(int k = 0; k < z; k++)
      {
        ind = i + x * (j + k * y);
        triplets.push_back(Eigen::Triplet<double>(ind, i + nx * (j + k * ny), 1));
        if (dir == 0) {
          triplets.push_back(Eigen::Triplet<double>(ind, (i + 1) + nx * (j + k * ny), -1));
        } else if (dir == 1) {
          triplets.push_back(Eigen::Triplet<double>(ind, i + nx * ((j + 1) + k * ny), -1));
        } else if (dir == 2) {
          triplets.push_back(Eigen::Triplet<double>(ind, i + nx * (j + (k + 1) * ny), -1));
        }
      }
    }
  }

  D.resize(x*y*z,nx*ny*nz);
  D.setZero();
  D.setFromTriplets(triplets.begin(), triplets.end());

}
