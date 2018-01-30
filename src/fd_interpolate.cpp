#include "fd_interpolate.h"
#include <cmath>
#include <vector>
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
  W.resize(P.rows(), nx * ny * nz);
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(8 * P.rows());
  for (int i = 0; i< P.rows(); i++){
    double x = std::floor((P(i, 0) - corner(0,0))/h); 
    double y = std::floor((P(i, 1) - corner(0,1))/h);
    double z = std::floor((P(i, 2) - corner(0,2))/h);
    double dx = (P(i, 0) - corner(0,0))/h - x;
    double dy = (P(i, 1) - corner(0,1))/h - y;
    double dz = (P(i, 2) - corner(0,2))/h - z;
    

    tripletList.push_back(Eigen::Triplet<double> (i, x + y *nx + z * nx*ny, (1-dx) * (1-dy) * (1-dz)));
    tripletList.push_back(Eigen::Triplet<double> (i, (x+1) + y *nx + z * nx*ny, dx * (1-dy) * (1-dz)));

    tripletList.push_back(Eigen::Triplet<double> (i, x + (y+1) *nx + z * nx*ny, (1-dx) * dy * (1-dz)));
    tripletList.push_back(Eigen::Triplet<double> (i, x + y * nx + (z+1) * nx*ny, (1-dx) * (1-dy) * dz));

    tripletList.push_back(Eigen::Triplet<double> (i, (x+1) + (y+1) *nx + z * nx*ny, dx * dy * (1-dz)));
    tripletList.push_back(Eigen::Triplet<double> (i, (x+1) + y *nx + (z+1) * nx*ny, dx  * (1-dy) * dz));

    tripletList.push_back(Eigen::Triplet<double> (i, x + (y+1) *nx + (z+1) * nx*ny, (1-dx) * dy * dz));
    tripletList.push_back(Eigen::Triplet<double> (i, (x+1) + (y+1) * nx + (z+1) * nx*ny, dx * dy * dz));

  }
  ////////////////////////////////////////////////////////////////////////////

  W.setFromTriplets(tripletList.begin(), tripletList.end());
}
