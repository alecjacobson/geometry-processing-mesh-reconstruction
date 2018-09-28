#include "fd_interpolate.h"
#include <iostream>
#include <math.h>

typedef Eigen::Triplet<double> T;

//lookup table
//http://paulbourke.net/miscellaneous/interpolation/
double get_coefficient(double dx, double dy, double dz, int corner)
{
  switch (corner)
  {
  case 0: //C000
    return (1 - dx) * (1 - dy) * (1 - dz);
  case 1: //C100
    return dx * (1 - dy) * (1 - dz);
  case 2: //C010
    return (1 - dx) * dy * (1 - dz);
  case 3: //C001
    return (1 - dx) * (1 - dy) * dz;
  case 4: //C101
    return dx * (1 - dy) * dz;
  case 5: //C011
    return (1 - dx) * dy * dz;
  case 6: //C110
    return dx * dy * (1 - dz);
  case 7: //C111
    return dx * dy * dz;
  }
}

void fd_interpolate(
    const int nx,
    const int ny,
    const int nz,
    const double h,
    const Eigen::RowVector3d &corner,
    const Eigen::MatrixXd &P,
    Eigen::SparseMatrix<double> &W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  //std::cout <<"corner is :"<< corner << std::endl;
  W.resize(P.rows(), nx * ny * nz);

  std::vector<T> tlist;
  tlist.reserve(8 * P.rows());

  for (int i = 0; i < P.rows(); i++)
  {
    double distance[] = {0.0, 0.0, 0.0};
    int grid[] = {0, 0, 0};
    double d[] = {0.0, 0.0, 0.0};

    for (int j = 0; j < 3; j++)
    {
      distance[j] = P(i, j) - corner[j];
      grid[j] = floor(distance[j] / h);
      d[j] = (distance[j] - grid[j] * h) / h;
    }

    tlist.push_back(T(i, grid[0] + grid[1] * nx + grid[2] * nx * ny,
                      get_coefficient(d[0], d[1], d[2], 0)));
    tlist.push_back(T(i, (grid[0] + 1) + grid[1] * nx + grid[2] * nx * ny,
                      get_coefficient(d[0], d[1], d[2], 1)));
    tlist.push_back(T(i, grid[0] + (grid[1] + 1) * nx + grid[2] * nx * ny,
                      get_coefficient(d[0], d[1], d[2], 2)));
    tlist.push_back(T(i, grid[0] + grid[1] * nx + (grid[2] + 1) * nx * ny,
                      get_coefficient(d[0], d[1], d[2], 3)));
    tlist.push_back(T(i, (grid[0] + 1) + grid[1] * nx + (grid[2] + 1) * nx * ny,
                      get_coefficient(d[0], d[1], d[2], 4)));
    tlist.push_back(T(i, grid[0] + (grid[1] + 1) * nx + (grid[2] + 1) * nx * ny,
                      get_coefficient(d[0], d[1], d[2], 5)));
    tlist.push_back(T(i, (grid[0] + 1) + grid[1] * nx + (grid[2] + 1) * nx * ny,
                      get_coefficient(d[0], d[1], d[2], 6)));
    tlist.push_back(T(i, (grid[0] + 1) + (grid[1] + 1) * nx + (grid[2] + 1) * nx * ny,
                      get_coefficient(d[0], d[1], d[2], 7)));

    /*
    W.coeffRef(i, grid[0] + grid[1] * nx + grid[2] * nx * ny) = get_coefficient(d[0],d[1],d[2], 0);
    W.coeffRef(i, (grid[0] + 1) + grid[1]* nx + grid[2] * nx * ny) = get_coefficient(d[0],d[1],d[2], 1);
    W.coeffRef(i, grid[0] + (grid[1] + 1) * nx + grid[2] * nx * ny) = get_coefficient(d[0],d[1],d[2], 2);
    W.coeffRef(i, grid[0] + grid[1] * nx + (grid[2] + 1) * nx * ny) = get_coefficient(d[0],d[1],d[2], 3);
    W.coeffRef(i, (grid[0] + 1) + grid[1] * nx + (grid[2] + 1) * nx * ny) = get_coefficient(d[0],d[1],d[2], 4);
    W.coeffRef(i, grid[0] + (grid[1] + 1) * nx + (grid[2] + 1) * nx * ny) = get_coefficient(d[0],d[1],d[2], 5);
    W.coeffRef(i, (grid[0] + 1) + grid[1] * nx + (grid[2] +1)* nx * ny) = get_coefficient(d[0],d[1],d[2], 6);
    W.coeffRef(i, (grid[0] + 1) + (grid[1] + 1) * nx + (grid[2] + 1) * nx * ny) = get_coefficient(d[0],d[1],d[2], 7);*/
  }

  //W.makeCompressed();
  W.setFromTriplets(tlist.begin(), tlist.end());
}
