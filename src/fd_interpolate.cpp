#include "fd_interpolate.h"
#include <math.h>
#include <vector>
#include <iostream>
typedef Eigen::Triplet<double> T;
void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  int n=P.rows();
  double x,y,z,dx,dy,dz;
  std::vector<T> tripletList;
  tripletList.reserve(n*8);
  W.resize(n,nx*ny*nz);
  for (int i = 0; i< n; i++){
    x=floor((P(i,0)-corner(0))/h); 
    y=floor((P(i,1)-corner(1))/h);
    z=floor((P(i, 2)-corner(2))/h);
    dx=(P(i,0)-corner(0))/h-x;
    dy=(P(i,1)-corner(1))/h-y;
    dz=(P(i,2)-corner(2))/h-z;

    tripletList.push_back(T(i, x + y *nx + z * nx*ny, (1-dx) * (1-dy) * (1-dz)));
    tripletList.push_back(T(i, (x+1) + y *nx + z * nx*ny, dx * (1-dy) * (1-dz)));

    tripletList.push_back(T(i, x + (y+1) *nx + z * nx*ny, (1-dx) * dy * (1-dz)));
    tripletList.push_back(T(i, x + y * nx + (z+1) * nx*ny, (1-dx) * (1-dy) * dz));

    tripletList.push_back(T(i, (x+1) + (y+1) *nx + z * nx*ny, dx * dy * (1-dz)));
    tripletList.push_back(T(i, (x+1) + y *nx + (z+1) * nx*ny, dx  * (1-dy) * dz));

    tripletList.push_back(T(i, x + (y+1) *nx + (z+1) * nx*ny, (1-dx) * dy * dz));
    tripletList.push_back(T(i, (x+1) + (y+1) * nx + (z+1) * nx*ny, dx * dy * dz));

  }
  W.setFromTriplets(tripletList.begin(),tripletList.end());
}
