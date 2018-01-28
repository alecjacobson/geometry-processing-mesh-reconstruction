#include "fd_interpolate.h"
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
  // for every p in P
    // Find the front-lower-left corner of the cube it belongs in
    // loop over the 8 points of the cube it belongs in
    // for every point in the cube, compute the area of the box with corners
      // assign the computed area to be the weight at W[p][7-idx of box]
  // h is 0.0217391
  std::vector< Eigen::Triplet<double> > tripletList;
  for(int p_idx = 0; p_idx < P.rows(); p_idx ++){
    // Find the front-lower-left corner of the cube it belongs in
    Eigen::RowVector3d point = Eigen::RowVector3d(P(p_idx,0), P(p_idx,1), P(p_idx,2));
    point -= corner;
    //std::cout << "Point: " << point << std::endl;
    Eigen::RowVector3d corner_of_p_cube;
    double x_cube = point(0) - std::fmod(point(0), h);
    double z_cube = point(1) - std::fmod(point(1), h);
    double y_cube = point(2) - std::fmod(point(2), h);
    corner_of_p_cube << x_cube , y_cube, z_cube;
    //std::cout << "Enclosing cube corner: " << corner_of_p_cube << std::endl;
    // get the x-offset, y-offset, z-offset into the cube
    float x_offset = point(0) - corner_of_p_cube(0);
    float y_offset = point(1) - corner_of_p_cube(1);
    float z_offset = point(2) - corner_of_p_cube(2);
    // find the i,j,k coordinates of the enclosing cube's corners
    int x_cube_idx = round(corner_of_p_cube(0) / h);
    int y_cube_idx = round(corner_of_p_cube(1) / h);
    int z_cube_idx = round(corner_of_p_cube(2) / h);
    // loop over the 8 points of the cube it belongs in
    for(int x_idx = 0; x_idx <= 1; x_idx += 1){
      for(int y_idx = 0; y_idx <= 1; y_idx += 1){
        for(int z_idx = 0; z_idx <= 1; z_idx += 1){
          // for every point in the cube, compute the area of the box with corners (x_cube, y_cube, z_cube) and P(P_idx)
          float x_area = (x_idx == 1) ? x_offset : 1 - x_offset;
          float y_area = (y_idx == 1) ? y_offset : 1 - y_offset;
          float z_area = (z_idx == 1) ? z_offset : 1 - z_offset;
          float area = x_area * y_area * z_area;
          // find the index of the current weight - the ijk coordinates of the vertex in the enclosing cube
          int i = x_cube_idx + x_idx;
          int j = y_cube_idx + y_idx;
          int k = z_cube_idx + z_idx;
          int weight_idx = i + nx*(j + k * ny);
          // assign the computed area to be the weight at W[p][idx_of_weight]
          if(weight_idx < W.cols()){
            tripletList.push_back(Eigen::Triplet<double>(p_idx, weight_idx, area));
          }
        }
      }
    }
  }
  W.setFromTriplets(tripletList.begin(), tripletList.end());
}
