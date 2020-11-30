#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
    typedef Eigen::Triplet<double> Tr;
    std::vector<Tr> trips;
    std::vector<int> n {nx, ny, nz};
    n[dir] -= 1;
    D.resize(n[0]*n[1]*n[2], nx*ny*nz);
    for (int i = 0; i < n[0]; i++) {
        for (int j = 0; j < n[1]; j++) {
            for (int k = 0; k < n[2]; k++) {
                trips.push_back(Tr(i + j*n[0] + k*n[0]*n[1], i + j*nx + k*nx*ny, -1));
                if (dir == 0) trips.push_back(Tr(i + j*n[0] + k*n[0]*n[1], (i+1) + j*nx + k*nx*ny, 1));
                else if (dir == 1) trips.push_back(Tr(i + j*n[0] + k*n[0]*n[1], i + (j+1)*nx + k*nx*ny, 1));
                else if (dir == 2) trips.push_back(Tr(i + j*n[0] + k*n[0]*n[1], i + j*nx + (k+1)*nx*ny, 1));
            }
        }
    }
    D.setFromTriplets(trips.begin(), trips.end());
}
