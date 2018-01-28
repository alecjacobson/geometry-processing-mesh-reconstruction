#include "fd_interpolate.h"

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
  int P_size = P.rows();
  for(int p_idx = 0; p_idx < P_size; p_idx ++){
    // Find the front-lower-left corner of the cube it belongs in
    Eigen::RowVector3d corner_of_p_cube;
    double x_cube = P(p_idx, 0) - std::fmod(P(p_idx, 0), h);
    double z_cube = P(p_idx, 1) - std::fmod(P(p_idx, 1), h);
    double y_cube = P(p_idx, 2) - std::fmod(P(p_idx, 2), h);
    corner_of_p_cube << x_cube , y_cube, z_cube;
    // get the x-offset, y-offset, z-offset into the cube
    float x_offset = P(p_idx, 0) - corner_of_p_cube(0);
    float y_offset = P(p_idx, 1) - corner_of_p_cube(1);
    float z_offset = P(p_idx, 2) - corner_of_p_cube(2);
    // loop over the 8 points of the cube it belongs in
    int weight_idx = 0;
    for(int x_idx = 0; x_idx <= 1; x_idx += 1){
      for(int y_idx = 0; y_idx <= 1; y_idx += 1){
        for(int z_idx = 0; z_idx <= 1; z_idx += 1){
          // for every point in the cube, compute the area of the box with corners (x_cube, y_cube, z_cube) and P(P_idx)
          float x_area = (x_idx == 1) ? x_offset : 1 - x_offset;
          float y_area = (y_idx == 1) ? y_offset : 1 - y_offset;
          float z_area = (z_idx == 1) ? z_offset : 1 - z_offset;
          float area = x_area * y_area * z_area;
          // find the index of the current weight
          // First, we must find x_cube, y_cube, z_cube - the coordinates of the vertex in the enclosing cube
          int x_cube = corner_of_p_cube(0) + (x_idx * h);
          int y_cube = corner_of_p_cube(1) + (y_idx * h);
          int z_cube = corner_of_p_cube(2) + (z_idx * h);
          // this is assuming that nx*ny*nz vertex positions are enumerated in the order they are visited via:
          // for x in nx; for y in ny; for z in nz; (x,y,z)
          //int weight_idx = (x_cube * ny*nz) + (y_cube*nz) + z_cube;
          // i + nx*(j + k * ny);
          int weight_idx = x_cube + nx*(y_cube + z_cube * ny);
          // assign the computed area to be the weight at W[p][idx_of_weight]
          W.insert(p_idx, weight_idx) = area;
          weight_idx += 1;
        }
      }
    }
  }
}
