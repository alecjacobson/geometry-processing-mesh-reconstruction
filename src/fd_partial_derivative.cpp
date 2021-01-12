#include "fd_partial_derivative.h"
#include <vector>
#include <iostream>

int find_col_idx_node(int idx_x, int idx_y, int idx_z, int nx, int ny);

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////

  // the number of rows
  int num_rows_D;

  // the actual matrix
  std::vector<Eigen::Triplet<double> > tripletList;

  // some intermediate variables
  int col_index_before, col_index_after, col_index;


  if (dir == 0) {
    num_rows_D = (nx-1)*ny*nz;
    // the upper bound for the number is 2 non-zero elements per row
    tripletList.reserve(num_rows_D*2);

    for (int i=0; i<nx-1; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz; k++) {
          // get the column index for the point and the next point
          col_index = find_col_idx_node(i, j, k, nx-1, ny);
          col_index_before = find_col_idx_node(i, j, k, nx, ny);
          col_index_after = find_col_idx_node(i+1, j, k, nx, ny);

          // assign the derivative
          tripletList.push_back(Eigen::Triplet<double>(col_index, col_index_before, -1.0));
          tripletList.push_back(Eigen::Triplet<double>(col_index, col_index_after, +1.0));
        }
      }
    }


  } else if (dir == 1) {
    num_rows_D = nx*(ny-1)*nz;
    // the upper bound for the number is 2 non-zero elements per row
    tripletList.reserve(num_rows_D*2);

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny-1; j++) {
        for (int k=0; k<nz; k++) {
          // get the column index for the point and the next point
          col_index = find_col_idx_node(i, j, k, nx, ny-1);
          col_index_before = find_col_idx_node(i, j, k, nx, ny);
          col_index_after = find_col_idx_node(i, j+1, k, nx, ny);

          // assign the derivative
          tripletList.push_back(Eigen::Triplet<double>(col_index, col_index_before, -1.0));
          tripletList.push_back(Eigen::Triplet<double>(col_index, col_index_after, +1.0));
        }
      }
    }

  } else {
    num_rows_D = nx*ny*(nz-1);
    // the upper bound for the number is 2 non-zero elements per row
    tripletList.reserve(num_rows_D*2);

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz-1; k++) {
          // get the column index for the point and the next point
          col_index_before = find_col_idx_node(i, j, k, nx, ny);
          col_index_after = find_col_idx_node(i, j, k+1, nx, ny);

          // assign the derivative
          tripletList.push_back(Eigen::Triplet<double>(col_index_before, col_index_before, -1.0));
          tripletList.push_back(Eigen::Triplet<double>(col_index_before, col_index_after, +1.0));
        }
      }
    }

  }
  

  
  D.resize(num_rows_D, nx*ny*nz);
  D.setZero();

  D.setFromTriplets(tripletList.begin(), tripletList.end());

  ////////////////////////////////////////////////////////////////////////////
}

// given the mesh node indices, returns the corresponding column index
int find_col_idx_node(int idx_x, int idx_y, int idx_z, int nx, int ny) {
  return idx_x + idx_y * nx + idx_z * ny * nx;
}