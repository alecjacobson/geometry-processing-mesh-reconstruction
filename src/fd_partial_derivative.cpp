#include "fd_partial_derivative.h"
#include <iostream>

// Idea Inspired by Eigen library manual of setFromTriplets() method
typedef Eigen::Triplet<double> T;

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  // Prepare new shape:
  int new_x = (dir == 0) ? nx - 1 : nx;
  int new_y = (dir == 1) ? ny - 1 : ny;
  int new_z = (dir == 2) ? nz - 1 : nz;
  
  int num_point = nx * ny * nz;
  int total = new_x * new_y * new_z;
  
  // Idea inspired by Eigen library manual of setFromTriplets() method
  std::vector<T> triplet_list; 

  // Each row contains two values (-1, 1):
  triplet_list.reserve(total*2);
  switch (dir) {
    case 0:
      for(int k = 0; k < new_z; k ++) {
        for(int j = 0; j < new_y; j ++) {
          for(int i = 0; i < new_x; i ++) {
            triplet_list.push_back(
              T(i + new_x*j + new_x*new_y*k,
                // skip previous nodes:
                i + nx*j + nx*ny*k, 
                -1/h));
            triplet_list.push_back(
              T(i + new_x*j + new_x*new_y*k,
                // move 1 column:
                i+1 + nx*j + nx*ny*k,
                1/h));
          }
        }
      }
      break;
    case 1:
      for(int k = 0; k < new_z; k ++) {
        for(int j = 0; j < new_y; j ++) {
          for(int i = 0; i < new_x; i ++) {
            triplet_list.push_back(
              T(i + new_x*j + new_x*new_y*k,
                // skip previous nodes:
                i + nx*j + nx*ny*k, 
                -1/h));
            triplet_list.push_back(
              T(i + new_x*j + new_x*new_y*k,
                // move 1 row:
                i + nx*(j+1) + nx*ny*k,
                1/h));
          }
        }
      }
      break;
    case 2:
      for(int k = 0; k < new_z; k ++) {
        for(int j = 0; j < new_y; j ++) {
          for(int i = 0; i < new_x; i ++) {
            triplet_list.push_back(
              T(i + new_x*j + new_x*new_y*k,
                // skip previous nodes:
                i + nx*j + nx*ny*k, 
                -1/h));
            triplet_list.push_back(
              T(i + new_x*j + new_x*new_y*k,
                // move 1 plane:
                i + nx*j + nx*ny*(k+1),
                1/h));
          }
        }
      }
      break;
  }

  D.resize(total, num_point);
  D.setFromTriplets(triplet_list.begin(), triplet_list.end());
  D.makeCompressed();
}
