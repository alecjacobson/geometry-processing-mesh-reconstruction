#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  switch(dir) {
  case 0: D.resize((nx-1)*ny*nz, nx*ny*nz); break;
  case 1: D.resize(nx*(ny-1)*nz, nx*ny*nz); break;
  case 2: D.resize(nx*ny*(nz-1), nx*ny*nz); break;
  }
  typedef Eigen::Triplet<double> T;
  std::vector<T> tps;
  tps.reserve(nx*ny*nz*2);
  for (int x = 0 ; x < nx ; ++x)
  for (int y = 0 ; y < ny ; ++y)
  for (int z = 0 ; z < nz ; ++z) {
    switch (dir) {
    case 0: if (x + 1 < nx) {
      tps.push_back(T(x*ny*nz+y*nz+z, x*ny*nz+y*nz+z, -1/h));
      tps.push_back(T(x*ny*nz+y*nz+z, (x+1)*ny*nz+y*nz+z, 1/h));
    } break;
    case 1: if (y + 1 < ny) {
      tps.push_back(T(x*(ny-1)*nz+y*nz+z, x*ny*nz+y*nz+z, -1/h));
      tps.push_back(T(x*(ny-1)*nz+y*nz+z, x*ny*nz+(y+1)*nz+z, 1/h));
    } break;
    case 2: if (z + 1 < nz) {
      tps.push_back(T(x*ny*(nz-1)+y*(nz-1)+z, x*ny*nz+y*nz+z, -1/h));
      tps.push_back(T(x*ny*(nz-1)+y*(nz-1)+z, x*ny*nz+y*nz+(z+1), 1/h));
    } break;
    }
  }
  D.setFromTriplets(tps.begin(), tps.end());
}
