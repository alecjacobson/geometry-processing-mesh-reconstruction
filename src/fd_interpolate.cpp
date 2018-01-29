#include "fd_interpolate.h"
#include <math.h>
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
  // number of input points
  const int n = P.rows();

  typedef Eigen::Triplet<double> tuple;
  std::vector<tuple> tupleList;
  tupleList.reserve(8*n);

  for (int m = 0; m < n; m++)
  {
    //Shift the volume to the origin
    double X = P(m, 0) - corner(0);
    double Y = P(m, 1) - corner(1);
    double Z = P(m, 2) - corner(2);

    //calculate grid coordinate
    int i = (int)trunc(X / h);
    int j = (int)trunc(Y / h);
    int k = (int)trunc(Z / h);

    //Normalize the single grid size to one
    double x = (X - i*h) / h;
    double y = (Y - j*h) / h;
    double z = (Z - k*h) / h;

    //Using the equations on http://paulbourke.net/miscellaneous/interpolation/
    //to calculate the weight of each vertices on the unit grid cube
    for(int r=0;r<=1;r++)
      for (int s = 0; s <= 1; s++)
        for (int t = 0; t <= 1; t++)
          tupleList.push_back(tuple(m, (i+r)+(j+s)*nx+(k+t)*nx*ny, (r*x+(1-r)*(1-x))*(s*y+(1-s)*(1-y))*(t*z+(1-t)*(1-z))));
  }
  W.setFromTriplets(tupleList.begin(), tupleList.end());
}
