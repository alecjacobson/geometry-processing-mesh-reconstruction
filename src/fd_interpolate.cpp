#include "fd_interpolate.h"
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
  W.resize(P.rows(), nx*ny*nz);
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(P.rows()*8);

  for(int i = 0; i < P.rows(); i += 1){

    int x = floor((P(i, 0) - corner(0))/h);
    int y = floor((P(i, 1) - corner(1))/h);
    int z = floor((P(i, 2) - corner(2))/h);

    double dx = (P(i, 0) - corner(0))/h - x;
    double dy = (P(i, 1) - corner(1))/h - y;
    double dz = (P(i, 2) - corner(2))/h - z;

    tripletList.push_back(T(i, x + nx*y + nx*ny*z, (1-dx)*(1-dy)*(1-dz)));
    tripletList.push_back(T(i, x + nx*y + nx*ny*(z+1), (1-dx)*(1-dy)*dz));
    tripletList.push_back(T(i, x + nx*(y+1) + nx*ny*z, (1-dx)*dy*(1-dz)));
    tripletList.push_back(T(i, x + nx*(y+1) + nx*ny*(z+1), (1-dx)*dy*dz));
    tripletList.push_back(T(i, (x+1) + nx*y + nx*ny*z, dx*(1-dy)*(1-dz)));
    tripletList.push_back(T(i, (x+1) + nx*y + nx*ny*(z+1), dx*(1-dy)*dz));
    tripletList.push_back(T(i, (x+1) + nx*(y+1) + nx*ny*z, dx*dy*(1-dz)));
    tripletList.push_back(T(i, (x+1) + nx*(y+1) + nx*ny*(z+1), dx*dy*dz));

  }
  W.setFromTriplets(tripletList.begin(), tripletList.end());
}
