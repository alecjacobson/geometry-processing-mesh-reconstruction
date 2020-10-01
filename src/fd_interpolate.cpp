#include "fd_interpolate.h"
#include <iostream>
#include <cmath>

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
  
  std::cout << "Building interpolation matrix..." << std::flush;

  // ****** DEBUGGING ******
  // std::cout << "Rows: " << P.rows() << ", Cols: " << P.cols() << std::endl;
  // std::cout << P(0,0) << ", " << P(0,1) << ", " << P(0,2) << std::endl;
  // std::cout << corner(0) << ", " << corner(1) << ", " << corner(2) << std::endl;
  // ***********************

  // Compute total number of grid points
  int ng = nx*ny*nz;

  // Extract total number of given points
  int np = P.rows();

  // Since the interpolation matrix multiplies into the list of grid points,
  // and interpolates them on to the given points, we know that the size of the
  // interpolation matrix must be np x ng. Also, since there are three components,
  // we extend the W matrix to concatenate each component on top of each other.

  W.resize(3*np, ng);
  W.reserve(3*np*4);

  // ****** DEBUGGING ******
  // std::cout << W.nonZeros() << std::endl;
  // ***********************

  // For each point in P, we need to find the nearest grid points, and then
  // use them to construct W.

  // Using the triplet insertion method based on Eigen's documentation
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;

  for (int ii = 0; ii < np; ii++)
  {
    // Compute distances from this point to the corner
    double dist_p_to_corner_x = std::abs(P(ii, 0) - corner(0));
    double dist_p_to_corner_y = std::abs(P(ii, 1) - corner(1));
    double dist_p_to_corner_z = std::abs(P(ii, 2) - corner(2));

    // Divide this distance by grid spacing - this is the number of cells between the point and the corner
    double cells_x = dist_p_to_corner_x/h;
    double cells_y = dist_p_to_corner_y/h;
    double cells_z = dist_p_to_corner_z/h;

    // By rounding these numbers up and down, we can figure out the nearest vertices
    double x_low = std::floor(cells_x);//*h + corner(0);
    double x_high = std::ceil(cells_x);//*h + corner(0);
    double y_low = std::floor(cells_y);//*h + corner(1);
    double y_high = std::ceil(cells_y);//*h + corner(1);
    double z_low = std::floor(cells_z);//*h + corner(2);
    double z_high = std::ceil(cells_z);//*h + corner(2);

    // Also figure out the *global* indices of the eight surrounding vertices
    int i000 = x_low + y_low*nx + z_low*nx*ny;
    int i100 = x_high + y_low*nx + z_low*nx*ny;
    int i010 = x_low + y_high*nx + z_low*nx*ny;
    int i110 = x_high + y_high*nx + z_low*nx*ny;
    int i001 = x_low + y_low*nx + z_high*nx*ny;
    int i101 = x_high + y_low*nx + z_high*nx*ny;
    int i011 = x_low + y_high*nx + z_high*nx*ny;
    int i111 = x_high + y_high*nx + z_high*nx*ny;

    // ****** DEBUGGING ******
    // std::cout << i000 << ", " <<  i100 << ", " << i010 << ", " << i110 << ", " << i001 << ", " << i101 << ", " << i011 << ", " << i111 << std::endl;  
    // ***********************

    // Now insert weights into appropriate spots of the interpolation matrix
    double w_x_low = 1 - (P(ii, 0) - x_low*h + corner(0))/(x_high*h + corner(0) - x_low*h + corner(0));
    double w_x_high = (P(ii, 0) - x_low*h + corner(0))/(x_high*h + corner(0) - x_low*h + corner(0));
    double w_y_low = 1 - (P(ii, 1) - y_low*h + corner(1))/(y_high*h + corner(1) - y_low*h + corner(1));
    double w_y_high = (P(ii, 1) - y_low*h + corner(1))/(y_high*h + corner(1) - y_low*h + corner(1));
    double w_z_low = 1 - (P(ii, 2) - z_low*h + corner(2))/(z_high*h + corner(2) - z_low*h + corner(2));
    double w_z_high = (P(ii, 2) - z_low*h + corner(2))/(z_high*h + corner(2) - z_low*h + corner(2));

    // x-weights
    tripletList.push_back(T(ii, i000, w_x_low));
    tripletList.push_back(T(ii, i010, w_x_low));
    tripletList.push_back(T(ii, i001, w_x_low));
    tripletList.push_back(T(ii, i011, w_x_low));

    tripletList.push_back(T(ii, i100, w_x_high));
    tripletList.push_back(T(ii, i110, w_x_high));
    tripletList.push_back(T(ii, i101, w_x_high));
    tripletList.push_back(T(ii, i111, w_x_high));

    // y-weights
    tripletList.push_back(T(ii + np, i000, w_y_low));
    tripletList.push_back(T(ii + np, i001, w_y_low));
    tripletList.push_back(T(ii + np, i101, w_y_low));
    tripletList.push_back(T(ii + np, i100, w_y_low));

    tripletList.push_back(T(ii + np, i010, w_y_high));
    tripletList.push_back(T(ii + np, i110, w_y_high));
    tripletList.push_back(T(ii + np, i011, w_y_high));
    tripletList.push_back(T(ii + np, i111, w_y_high));

    // z-weights
    tripletList.push_back(T(ii + 2*np, i000, w_z_low));
    tripletList.push_back(T(ii + 2*np, i100, w_z_low));
    tripletList.push_back(T(ii + 2*np, i010, w_z_low));
    tripletList.push_back(T(ii + 2*np, i110, w_z_low));

    tripletList.push_back(T(ii + 2*np, i001, w_z_high));
    tripletList.push_back(T(ii + 2*np, i101, w_z_high));
    tripletList.push_back(T(ii + 2*np, i011, w_z_high));
    tripletList.push_back(T(ii + 2*np, i111, w_z_high));

  
    // ****** DEBUGGING ******
    // std::cout << x_low << ", " << P(ii, 0) << ", " << x_high << std::endl;
    // std::cout << y_low << ", " << P(ii, 1) << ", " << y_high << std::endl;
    // std::cout << z_low << ", " << P(ii, 2) << ", " << z_high << std::endl;

    // std::cout << (P(ii, 0) - x_low)/(x_high - x_low) << std::endl;
    // std::cout << (P(ii, 1) - y_low)/(y_high - y_low) << std::endl;
    // std::cout << (P(ii, 2) - z_low)/(z_high - z_low) << std::endl;
    // ***********************

  }

  W.setFromTriplets(tripletList.begin(), tripletList.end());
  W.finalize();

  std::cout << "done." << std::endl;

  ////////////////////////////////////////////////////////////////////////////
}


