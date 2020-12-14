#include "fd_interpolate.h"
#include <iostream>

using namespace std;
using namespace Eigen;

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  vector<Triplet<double>> triplets;
  triplets.reserve(P.rows() * 8);

  for (int i = 0; i < P.rows(); i++) {
    RowVector3d diff = P.row(i) - corner;

    double xr = diff(0) * 1.0 / h;
    double yr = diff(1) * 1.0 / h;
    double zr = diff(2) * 1.0 / h;

    int xi = floor(xr);
    double xd = xr - xi;

    int yi = floor(yr);
    double yd = yr - yi;

    int zi = floor(zr);
    double zd = zr - zi;

    // xi, yi, zi
    triplets.push_back(Triplet<double>(i, nx * ny * zi + nx * yi + xi, (1 - xd) * (1 - yd) * (1 - zd)));
    // xi + 1, yi, zi
    triplets.push_back(Triplet<double>(i, nx * ny * zi + nx * yi + (xi + 1), xd * (1 - yd) * (1 - zd)));
    // xi, yi + 1, zi
    triplets.push_back(Triplet<double>(i, nx * ny * zi + nx * (yi + 1) + xi, (1 - xd) * yd * (1 - zd)));
    // xi, yi, zi + 1
    triplets.push_back(Triplet<double>(i, nx * ny * (zi + 1) + nx * yi + xi, (1 - xd) * (1 - yd) * zd));
    // xi + 1, yi + 1, zi
    triplets.push_back(Triplet<double>(i, nx * ny * zi + nx * (yi + 1) + (xi + 1), xd * yd * (1 - zd)));
    // xi + 1, yi, zi + 1
    triplets.push_back(Triplet<double>(i, nx * ny * (zi + 1) + nx * yi + (xi + 1), xd * (1 - yd) * zd));
    // xi, yi + 1, zi + 1
    triplets.push_back(Triplet<double>(i, nx * ny * (zi + 1) + nx * (yi + 1) + xi, (1 - xd) * yd * zd));
    // xi + 1, yi + 1, zi + 1
    triplets.push_back(Triplet<double>(i, nx * ny * (zi + 1) + nx * (yi + 1) + (xi + 1), xd * yd * zd));
  }

  W.resize(P.rows(), nx*ny*nz);
  W.setFromTriplets(triplets.begin(), triplets.end());
}
