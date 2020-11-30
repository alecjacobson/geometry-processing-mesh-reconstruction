#include "fd_partial_derivative.h"
#include <vector>

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  int nx_p = nx - (dir == 0);
  int ny_p = ny - (dir == 1);
  int nz_p = nz - (dir == 2);

  D.resize(nx_p*ny_p*nz_p, nx*ny*nz);
  
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  
  for (int i = 0; i < nx_p; i++) {
    for (int j = 0; j < ny_p; j++) {
      for (int k = 0; k < nz_p; k++) {
        int row = i + nx_p*(j + k * ny_p);
        int left_ind = i + nx*(j + k * ny);
        int right_ind = i + (dir == 0) + nx*(j + (dir == 1) + (k + (dir == 2))* ny);
        tripletList.push_back(T(row,left_ind,-1/h));
        tripletList.push_back(T(row,right_ind,1/h));
      }
    }
  }
  
  D.setFromTriplets(tripletList.begin(), tripletList.end());
}
