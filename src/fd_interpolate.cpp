#include "fd_interpolate.h"

#include <iostream>
#include <vector>

#include "igl/floor.h"

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  // Get the grid corner indices for each input point
  Eigen::MatrixXd P_rel_grid = (P.rowwise() - corner) / h;

  Eigen::MatrixXi P_corners_ijk;
  igl::floor(P_rel_grid, P_corners_ijk);

  // These are the first three weights
  Eigen::MatrixXd P_rel_corner = P_rel_grid - P_corners_ijk.cast<double>();

  auto ind = [&nx, &ny](int i, int j, int k) {
    return i + nx * (j + k * ny);
  };

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(8 * P.rows());
  for(int l = 0; l < P.rows(); l++) 
  {
    int i = P_corners_ijk(l, 0);
    int j = P_corners_ijk(l, 1);
    int k = P_corners_ijk(l, 2);

    double wx = P_rel_corner(l, 0);
    double wy = P_rel_corner(l, 1);
    double wz = P_rel_corner(l, 2);

    tripletList.push_back(T(l, ind(i+0, j+0, k+0), (1. - wx) * (1. - wy) * (1. - wz)));
    tripletList.push_back(T(l, ind(i+0, j+0, k+1), (1. - wx) * (1. - wy) * wz));
    tripletList.push_back(T(l, ind(i+0, j+1, k+0), (1. - wx) * wy        * (1. - wz)));
    tripletList.push_back(T(l, ind(i+0, j+1, k+1), (1. - wx) * wy        * wz));
    tripletList.push_back(T(l, ind(i+1, j+0, k+0), wx        * (1. - wy) * (1. - wz)));
    tripletList.push_back(T(l, ind(i+1, j+0, k+1), wx        * (1. - wy) * wz));
    tripletList.push_back(T(l, ind(i+1, j+1, k+0), wx        * wy        * (1. - wz)));
    tripletList.push_back(T(l, ind(i+1, j+1, k+1), wx        * wy        * wz));
  }

  W.resize(P.rows(), nx*ny*nz);
  W.setFromTriplets(tripletList.begin(), tripletList.end());
}
