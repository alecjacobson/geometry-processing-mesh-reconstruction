#include "fd_interpolate.h"
#include <cmath>
#include <vector>
#include <iostream>

int compute_col_idx_node(int idx_x, int idx_y, int idx_z, int nx, int ny);

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

  // Assumption: The Point list P is of dimensions: integer X 3 (That is, coordinates represented by columns).

  // For each of the points, need to find the 1-step grid range in x, y, z to localize the corresponding cube
  // that is required for trilinear interpolation.

  // requried variables that localize the small unit cube in the grid for a given point 
  double grid_cube_x_0;
  double grid_cube_x_1;
  double grid_cube_y_0;
  double grid_cube_y_1;
  double grid_cube_z_0;
  double grid_cube_z_1;

  // for indexing the grid nodes
  int grid_cube_x_0_idx;
  int grid_cube_x_1_idx;
  int grid_cube_y_0_idx;
  int grid_cube_y_1_idx;
  int grid_cube_z_0_idx;
  int grid_cube_z_1_idx;

  // define intermediate variables
  double x_d, y_d, z_d;
  double c_000, c_001, c_100, c_010, c_110, c_101, c_011, c_111;


  // the actual matrix
  std::vector<Eigen::Triplet<double> > tripletList;

  // W would be of dimensions num_points X (nx*ny*nz)
  // But each row would have utmost 8 non-zero entries
  // Therefore, reserve the Triplet size to num_points * 8
  tripletList.reserve(P.rows()*8);



  for (int i=0; i<P.rows(); i++) {
    // for each point given, we need the corresponding small cube of the grid to be identified 
    // as per the question, the corner of the bottom-left-front-most node of the grid is given

    // grid_cube_x_0, grid_cube_x_1 represent the two x-coordinates between which the x-coordinate of the point exists
    grid_cube_x_0 = corner(0) + double(floor((P(i,0)-corner(0))/h)) * h;
    grid_cube_x_1 = grid_cube_x_0 + h;

    // similarly for the other coordinates
    grid_cube_y_0 = corner(1) + double(floor((P(i,1)-corner(1))/h)) * h;
    grid_cube_y_1 = grid_cube_y_0 + h;

    grid_cube_z_0 = corner(2) + double(floor((P(i,2)-corner(2))/h)) * h;
    grid_cube_z_1 = grid_cube_z_0 + h;

    // now we have localized the cube to which the point p(i) belong to
    // let's assign integers to the nodes in the mesh so that it is easy to vectorize as X
    // this helps in knowing where to fill in the W matrix
    // Taking Corner index as (0,0,0), we assign indices to all nodes assuming that the indices are all positive (since we need to access it as a matrix index)
    grid_cube_y_0_idx = int((grid_cube_y_0 - corner(1))/h);
    grid_cube_y_1_idx = grid_cube_y_0_idx + 1;

    grid_cube_z_0_idx = int((grid_cube_z_0 - corner(2))/h);
    grid_cube_z_1_idx = grid_cube_z_0_idx + 1;

    grid_cube_x_0_idx = int((grid_cube_x_0 - corner(0))/h);
    grid_cube_x_1_idx = grid_cube_x_0_idx + 1;



    // interpolation matrix W's weights are going to based on the distances between the above identified small cube
    // and the given point - so, W is of dimensions: num_points X (nx * ny * nz)
    // and given (i, j, k) integer mesh-coordinates of a node, we have it corresponding 
    // to the column index: i + j * n_x + k * n_y * n_x - in W 
    // compute the interpolation weights wrt to the cube identified for the point
    
    // define intermediate variables for computing interpolation
    x_d = (P(i,0) - grid_cube_x_0)/h;
    y_d = (P(i,1) - grid_cube_y_0)/h;
    z_d = (P(i,2) - grid_cube_z_0)/h;

    // compute the weights
    c_000 = (1.0 - x_d) * (1.0 - y_d) * (1.0 - z_d);
    c_100 = (x_d) * (1.0 - y_d) * (1.0 - z_d);
    c_010 = (1.0 - x_d) * (y_d) * (1.0 - z_d);
    c_110 = (x_d) * (y_d) * (1.0 - z_d);
    c_001 = (1.0 - x_d) * (1.0 - y_d) * (z_d);
    c_101 = (x_d) * (1.0 - y_d) * (z_d);
    c_011 = (1.0 - x_d) * (y_d) * (z_d);
    c_111 = (x_d) * (y_d) * (z_d);



    // insert the 8 weights at appropriate columns
    tripletList.push_back(Eigen::Triplet<double>(i, compute_col_idx_node(grid_cube_x_0_idx, grid_cube_y_0_idx, grid_cube_z_0_idx, nx, ny), c_000));
    tripletList.push_back(Eigen::Triplet<double>(i, compute_col_idx_node(grid_cube_x_1_idx, grid_cube_y_0_idx, grid_cube_z_0_idx, nx, ny), c_100));
    tripletList.push_back(Eigen::Triplet<double>(i, compute_col_idx_node(grid_cube_x_0_idx, grid_cube_y_1_idx, grid_cube_z_0_idx, nx, ny), c_010));
    tripletList.push_back(Eigen::Triplet<double>(i, compute_col_idx_node(grid_cube_x_1_idx, grid_cube_y_1_idx, grid_cube_z_0_idx, nx, ny), c_110));
    tripletList.push_back(Eigen::Triplet<double>(i, compute_col_idx_node(grid_cube_x_0_idx, grid_cube_y_0_idx, grid_cube_z_1_idx, nx, ny), c_001));
    tripletList.push_back(Eigen::Triplet<double>(i, compute_col_idx_node(grid_cube_x_1_idx, grid_cube_y_0_idx, grid_cube_z_1_idx, nx, ny), c_101));
    tripletList.push_back(Eigen::Triplet<double>(i, compute_col_idx_node(grid_cube_x_0_idx, grid_cube_y_1_idx, grid_cube_z_1_idx, nx, ny), c_011));
    tripletList.push_back(Eigen::Triplet<double>(i, compute_col_idx_node(grid_cube_x_1_idx, grid_cube_y_1_idx, grid_cube_z_1_idx, nx, ny), c_111));


    // iterate
  }



W.resize(P.rows(), nx*ny*nz);
W.setZero();

// set the values to the sparse matrices
W.setFromTriplets(tripletList.begin(), tripletList.end());

  ////////////////////////////////////////////////////////////////////////////
}

// given the mesh node indices, returns the column index in W that corresponds to its interpolation weight contribution
int compute_col_idx_node(int idx_x, int idx_y, int idx_z, int nx, int ny) {
  int ind = idx_x + idx_y * nx + idx_z * ny * nx;
  return ind;
}
