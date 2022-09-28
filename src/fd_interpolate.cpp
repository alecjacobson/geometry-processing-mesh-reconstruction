#include "fd_interpolate.h"
#include <iostream>

// Idea Inspired by Eigen library manual of setFromTriplets() method
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
  // Construct a matrix(P.rows, P.cols) with all rows equal to corner
  Eigen::MatrixXd corner_matrix(P.rows(), P.cols());
  corner_matrix << Eigen::MatrixXd::Ones(P.rows(), P.cols());

  corner_matrix.col(0) << corner_matrix.col(0) * corner.col(0);
  corner_matrix.col(1) << corner_matrix.col(1) * corner.col(1);
  corner_matrix.col(2) << corner_matrix.col(2) * corner.col(2);

  // Compute relative position in the grid
  Eigen::MatrixXd position(P.rows(), P.cols());
  position << (P - corner_matrix) / h;

  Eigen::MatrixXd grid(P.rows(), P.cols());
  grid << (position.cast <int> ()).cast <double> ();

  Eigen::MatrixXd abso(P.rows(), P.cols());
  abso << position - grid;


  // Idea inspired by Eigen library manual of setFromTriplets() method
  std::vector<T> triplet_list; 
  triplet_list.reserve(P.rows() * 8);

  for (int i = 0; i < P.rows(); i ++) {
    // x, y, z
    triplet_list.push_back(
      T(i, grid(i, 0) + grid(i, 1)*nx + grid(i, 2)*nx*ny, (1-abso(i, 0))*(1-abso(i, 1))*(1-abso(i, 2)) )
    );
    // x+1, y, z
    triplet_list.push_back(
      T(i, grid(i, 0)+1 + grid(i, 1)*nx + grid(i, 2)*nx*ny, abso(i, 0)*(1-abso(i, 1))*(1-abso(i, 2)))
    );
    // x, y+1, z
    triplet_list.push_back(
      T(i, grid(i, 0) + (grid(i, 1)+1)*nx + grid(i, 2)*nx*ny,  (1-abso(i, 0))*abso(i,1)*(1-abso(i, 2)))
    );
    // x+1, y+1, z
    triplet_list.push_back(
      T(i, grid(i, 0)+1 + (grid(i, 1)+1)*nx + grid(i, 2)*nx*ny, abso(i, 0)*abso(i, 1)*(1-abso(i, 2)))
    );
    // x, y, z+1
    triplet_list.push_back(
      T(i, grid(i, 0) + grid(i, 1)*nx + (grid(i, 2)+1)*nx*ny, (1-abso(i, 0))*(1-abso(i, 1))*abso(i, 2) )
    );
    // x+1, y, z+1
    triplet_list.push_back(
      T(i, grid(i, 0)+1 + grid(i, 1)*nx + (grid(i, 2)+1)*nx*ny, abso(i, 0)*(1-abso(i, 1))*abso(i, 2))
    );
    // x, y+1, z+1
    triplet_list.push_back(
      T(i, grid(i, 0) + (grid(i, 1)+1)*nx + (grid(i, 2)+1)*nx*ny, (1-abso(i, 0))*abso(i, 1)*abso(i, 2))
    );
    // x+1, y+1, z+1
    triplet_list.push_back(
      T(i, grid(i, 0)+1 + (grid(i, 1)+1)*nx + (grid(i, 2)+1)*nx*ny, abso(i, 0)*abso(i, 1)*abso(i, 2))
    );
    
  }
  W.resize(P.rows(), nx*ny*nz);
  W.setFromTriplets(triplet_list.begin(), triplet_list.end());
  W.makeCompressed();
}
