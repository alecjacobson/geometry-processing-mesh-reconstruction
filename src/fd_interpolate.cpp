#include "fd_interpolate.h"
#include <cmath>
#include <iostream>

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  W.resize(P.rows(),nx * ny * nz);
  typedef Eigen::Triplet<double> T;
  std::vector< T > tripletList;
  tripletList.reserve(8*P.rows());

  for (int i = 0; i < P.rows(); i++) {
    // distance
    double dx = P(i,0) - corner[0];
    double dy = P(i,1) - corner[1];
    double dz = P(i,2) - corner[2];

    // position of grids
    int gx = floor(dx / h);
    int gy = floor(dy / h);
    int gz = floor(dz / h);

    // relative distance
    double rx = (dx - gx*h) / h;
    double ry = (dy - gy*h) / h ;
    double rz = (dz - gz*h) / h;

    // eight point of one cube
    tripletList.push_back(T(i, gx+nx*(gy+gz*ny), (1-rx)*(1-ry)*(1-rz)));
    tripletList.push_back(T(i, gx+1+nx*(gy+gz*ny), rx*(1-ry)*(1-rz)));
    tripletList.push_back(T(i, gx+nx*(gy+1+gz*ny), (1-rx)*ry*(1-rz)));

    tripletList.push_back(T(i, gx+nx*(gy+(gz+1)*ny), (1-rx)*(1-ry)*rz));
    tripletList.push_back(T(i, gx+1+nx*(gy+1+gz*ny), rx*ry*(1-rz)));
    tripletList.push_back(T(i, gx+1+nx*(gy+(gz+1)*ny), rx*(1-ry)*rz));

    tripletList.push_back(T(i, gx+nx*(gy+1+(gz+1)*ny), (1-rx)*ry*rz));
    tripletList.push_back(T(i, gx+1+nx*(gy+1+(gz+1)*ny), rx*ry*rz));

  }
  W.setFromTriplets(tripletList.begin(), tripletList.end());
}
