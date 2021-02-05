#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
    typedef Eigen::Triplet<double> Triplets;

    auto to_idx = [=](int i, int j, int k) {
        return i + nx * (j + ny * k);
    };

    auto to_sidx = [=](int i, int j, int k) {
        return i + (nx - (dir==0 ? 1 : 0)) * (j + (ny - (dir==1 ? 1 : 0)) * k);
    };

    int next_idx = 0;
    std::vector<Triplets> triplets;
    triplets.reserve(nx*ny*nz*2);

    for (int i = 0; i < nx - (dir==0 ? 1 : 0); ++i) {
        for (int j = 0; j < ny - (dir==1 ? 1 : 0); ++j) {
            for (int k = 0; k < nz - (dir==2 ? 1 : 0); ++k) {
                if (dir == 0)
                    next_idx = to_idx(i+1, j, k);
                if (dir == 1)
                    next_idx = to_idx(i, j+1, k);
                if (dir == 2)
                    next_idx = to_idx(i, j, k+1);
                triplets.emplace_back(to_sidx(i, j, k), to_idx(i, j, k), -1);
                triplets.emplace_back(to_sidx(i, j, k), next_idx, 1);
            }
        }
    }

    D.setFromTriplets(triplets.begin(), triplets.end());
}
