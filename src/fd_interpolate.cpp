#include "fd_interpolate.h"

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  W.resize(P.rows(), (nx*ny*nz));
  typedef Eigen::Triplet<double> T;
  std::vector<T> tps;
  tps.reserve(P.rows()*8);
  for ( int i = 0 ; i < P.rows() ; ++i ) {
    Eigen::RowVector3d p(P(i, 0), P(i, 1), P(i, 2));
    Eigen::RowVector3d g = (p - corner) / h;
    int ix = int(g[0]);
    int iy = int(g[1]);
    int iz = int(g[2]);
    double cx = g[0] - ix;
    double cy = g[1] - iy;
    double cz = g[2] - iz;
    tps.push_back(T(i, ix+iy*nx+iz*nx*ny, (1-cx)*(1-cy)*(1-cz)));
    tps.push_back(T(i, ix+iy*nx+(iz+1)*nx*ny, (1-cx)*(1-cy)*cz));
    tps.push_back(T(i, ix+(iy+1)*nx+iz*nx*ny, (1-cx)*cy*(1-cz)));
    tps.push_back(T(i, ix+(iy+1)*nx+(iz+1)*nx*ny, (1-cx)*cy*cz));
    tps.push_back(T(i, (ix+1)+iy*nx+iz*nx*ny, cx*(1-cy)*(1-cz)));
    tps.push_back(T(i, (ix+1)+iy*nx+(iz+1)*nx*ny, cx*(1-cy)*cz));
    tps.push_back(T(i, (ix+1)+(iy+1)*nx+iz*nx*ny, cx*cy*(1-cz)));
    tps.push_back(T(i, (ix+1)+(iy+1)*nx+(iz+1)*nx*ny, cx*cy*cz));

  }
  W.setFromTriplets(tps.begin(), tps.end());
}
