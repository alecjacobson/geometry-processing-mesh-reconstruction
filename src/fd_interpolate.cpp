#include "fd_interpolate.h"
#include "fd_get_ind.h"
#include <iostream>

using namespace std;

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  int num_rows_P = P.rows();

  std::vector<Eigen::Triplet<double>> triplets;
  // We have 8 non-zero values per row (namely, the cube surrounding a point)  
  triplets.reserve(8 * num_rows_P);
  for(int i = 0; i < num_rows_P; i++)
  {
    int x_pos, y_pos, z_pos;
    int c000, c001, c010, c011, c100, c101, c110, c111;
    double x_d, y_d, z_d;

    // Determine the point's position in the grid space
    Eigen::RowVector3d position = P.row(i);
    Eigen::RowVector3d grid_position = (position - corner) / h;

    x_pos = grid_position[0];
    y_pos = grid_position[1];
    z_pos = grid_position[2];

    // Indices corresponding to the cube surrounding this point
    c000 = fd_get_ind(nx, ny, nz, x_pos, y_pos, z_pos);
    c001 = fd_get_ind(nx, ny, nz, x_pos, y_pos, z_pos + 1);
    c010 = fd_get_ind(nx, ny, nz, x_pos, y_pos + 1, z_pos);
    c011 = fd_get_ind(nx, ny, nz, x_pos, y_pos + 1, z_pos + 1);
    c100 = fd_get_ind(nx, ny, nz, x_pos + 1, y_pos, z_pos);
    c101 = fd_get_ind(nx, ny, nz, x_pos + 1, y_pos, z_pos + 1);
    c110 = fd_get_ind(nx, ny, nz, x_pos + 1, y_pos + 1, z_pos);
    c111 = fd_get_ind(nx, ny, nz, x_pos + 1, y_pos + 1, z_pos + 1);

    // Differences
    x_d = grid_position[0] - x_pos;
    y_d = grid_position[1] - y_pos;
    z_d = grid_position[2] - z_pos;

    triplets.push_back(Eigen::Triplet<double>(i, c000, (1.0 - x_d)*(1.0 - y_d)*(1.0 - z_d)));
    triplets.push_back(Eigen::Triplet<double>(i, c001, (1.0 - x_d)*(1.0 - y_d)*(z_d)));
    triplets.push_back(Eigen::Triplet<double>(i, c010, (1.0 - x_d)*(y_d)*(1.0 - z_d)));
    triplets.push_back(Eigen::Triplet<double>(i, c011, (1.0 - x_d)*(y_d)*(z_d)));
    triplets.push_back(Eigen::Triplet<double>(i, c100, (x_d)*(1.0 - y_d)*(1.0 - z_d)));
    triplets.push_back(Eigen::Triplet<double>(i, c101, (x_d)*(1.0 - y_d)*(z_d)));
    triplets.push_back(Eigen::Triplet<double>(i, c110, (x_d)*(y_d)*(1.0 - z_d)));
    triplets.push_back(Eigen::Triplet<double>(i, c111, (x_d)*(y_d)*(z_d)));
  }
  W.setFromTriplets(triplets.begin(), triplets.end());
}
