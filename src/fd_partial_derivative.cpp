#include "fd_partial_derivative.h"

using namespace Eigen;
using namespace std;

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  int nxd = nx;
  int nyd = ny;
  int nzd = nz;

  if (dir == 0) {
    nxd--;
  }
  else if (dir == 1) {
    nyd--;
  }
  else if (dir == 2) {
    nzd--;
  }

  vector<Triplet<double>> triplets;
  triplets.reserve(D.rows()*2);

  for (int i = 0; i < nxd; i++) {
    for (int j = 0; j < nyd; j++) {
      for (int k = 0; k < nzd; k++) {
        // same for all the three cases
        int cur = nxd * nyd * k + nxd * j + i;
        triplets.push_back(Triplet<double>(cur, nx * ny * k + nx * j + i, - 1.0));

        if (dir == 0) {
          triplets.push_back(Triplet<double>(cur, nx * ny * k + nx * j + (i + 1), 1.0));
        }
        else if (dir == 1) {
          triplets.push_back(Triplet<double>(cur, nx * ny * k + nx * (j + 1) + i, 1.0));
        }
        else if (dir == 2) {
          triplets.push_back(Triplet<double>(cur, nx * ny * (k + 1) + nx * j + i, 1.0));
        }
      }
    }
  }

  D.resize(nxd*nyd*nzd, nx*ny*nz);
  D.setFromTriplets(triplets.begin(), triplets.end());
}
