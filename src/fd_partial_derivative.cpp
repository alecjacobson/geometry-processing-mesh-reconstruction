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
  ////////////////////////////////////////////////////////////////////////////

  // staggered dimensions
  int nx_ = nx;
  int ny_ = ny;
  int nz_ = nz;
  if (dir == 0)
    nx_ -= 1;
  else if (dir == 1)
    ny_ -= 1;
  else if (dir == 2)
    nz_ -= 1;
  else
    throw std::runtime_error("[fd_partial_derivative] dir must be 0, 1, or 2.");

  int m = nx_*ny_*nz_;

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(2*nx_*ny_*nz_);

  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      for (int k = 0; k < nz_; ++k) {
        int stag_id = i + j*nx_ + k*ny_*nx_;
        int l_ = i + j*nx + k*ny*nx;
        int l;
        if (dir == 0)
          l = i + 1 + j*nx + k*ny*nx;
        else if (dir == 1)
          l = i + (j + 1)*nx + k*ny*nx;
        else if (dir == 2)
          l = i + j*nx + (k + 1)*ny*nx;

        tripletList.push_back(T(stag_id, l, 1));
        tripletList.push_back(T(stag_id, l_, -1));
      } // end loop k
    } // end loop j
  } // end loop i

  // create sparse matrix
  D = Eigen::SparseMatrix<double>(m, nx*ny*nz);
  D.setFromTriplets(tripletList.begin(), tripletList.end());
}
