#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  int is_x = dir == 0;
  int is_y = dir == 1;
  int is_z = dir == 2;

  int nxd = nx - 1 * is_x;
  int nyd = ny - 1 * is_y;
  int nzd = nz - 1 * is_z;

  auto ind = [&nx, &ny](int i, int j, int k) {
    return i + nx * (j + k * ny);
  };

  auto ind_d = [&nxd, &nyd](int i, int j, int k) {
    return i + nxd * (j + k * nyd);
  };

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(nxd*nyd*nzd);

  for(int i = 0; i < nxd; i++) 
  {
    for(int j = 0; j < nyd; j++)
    {
      for(int k = 0; k < nzd; k++)
      {
        int ip = i + 1 * is_x;
        int jp = j + 1 * is_y;
        int kp = k + 1 * is_z;

        int row = ind_d(i, j, k);
        tripletList.push_back(T(row, ind(i,j,k),     1./h));
        tripletList.push_back(T(row, ind(ip,jp,kp), -1./h));
      }
    }
  }

  D.resize(nxd*nyd*nzd, nx*ny*nz);
  D.setFromTriplets(tripletList.begin(), tripletList.end());
}
