#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  G.resize((nx-1)*ny*nz+nx*(ny-1)*nz+nx*ny*(nz-1), nx*ny*nz);
  int s0 = 0;
  int s1 = (nx-1)*ny*nz;
  int s2 = s1+nx*(ny-1)*nz;
  typedef Eigen::Triplet<double> T;
  std::vector<T> tps;
  tps.reserve(nx*ny*nz*6);
  for (int x = 0 ; x < nx ; ++x)
  for (int y = 0 ; y < ny ; ++y)
  for (int z = 0 ; z < nz ; ++z) {
    if (x + 1 < nx) {
      tps.push_back(T(s0+x+y*(nx-1)+z*(nx-1)*ny, x+y*nx+z*nx*ny, -1/h));
      tps.push_back(T(s0+x+y*(nx-1)+z*(nx-1)*ny, (x+1)+y*nx+z*nx*ny, 1/h));
    }
    if (y + 1 < ny) {
      tps.push_back(T(s1+x+y*nx+z*nx*(ny-1), x+y*nx+z*nx*ny, -1/h));
      tps.push_back(T(s1+x+y*nx+z*nx*(ny-1), x+(y+1)*nx+z*nx*ny, 1/h));
    }
    if (z + 1 < nz) {
      tps.push_back(T(s2+x+y*nx+z*nx*ny, x+y*nx+z*nx*ny, -1/h));
      tps.push_back(T(s2+x+y*nx+z*nx*ny, x+y*nx+(z+1)*nx*ny, 1/h));
    }
  }
  G.setFromTriplets(tps.begin(), tps.end());
}
