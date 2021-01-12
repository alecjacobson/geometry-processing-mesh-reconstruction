#include "fd_interpolate.h"
#include <iostream>
#include <cmath>
#include <string>

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
  int i_int, j_int, k_int;
  double i, j, k;
  double w_x, w_y, w_z;
  // P -= corner;
  // Eigen::MatrixXd temp = ((P).rowwise() - corner) / h;
  // Eigen::MatrixXi P_int = temp.cast <int> ();
  // Eigen::MatrixXd P_frac = temp - P_int.cast<double> ();

  std::vector<Eigen::Triplet<double> > tripletList;

  for (int row=0; row<P.rows(); row++) {

    i = (P(row,0) - corner(0)) / h;
    j = (P(row,1) - corner(1)) / h;
    k = (P(row,2) - corner(2)) / h;

    i_int = std::floor(i);
    j_int = std::floor(j);
    k_int = std::floor(k);

    w_x = i - i_int;
    w_y = j - j_int;
    w_z = k - k_int;

    // const auto col = i_int + nx*(j_int + k_int * ny);
    // W.insert(row, col) = (w_x*w_y*w_z);
    // W.insert(row, col+1) = ((1-w_x)*w_y*w_z);
    // W.insert(row, col+nx) = (w_x*(1-w_y)*w_z);
    // W.insert(row, col+nx*ny) = (w_x*w_y*(1-w_z));
    // W.insert(row, col+1+nx) = ((1-w_x)*(1-w_y)*w_z);
    // W.insert(row, col+1+nx*ny) = ((1-w_x)*w_y*(1-w_z));
    // W.insert(row, col+nx+nx*ny) = (w_x*(1-w_y)*(1-w_z));
    // W.insert(row, col+1+nx+nx*ny) = ((1-w_x)*(1-w_y)*(1-w_z));

    const auto col = i_int + nx*(j_int + k_int * ny);

    // tripletList.push_back(Eigen::Triplet<double>(row, col, P_frac(row,0) * P_frac(row,1) * P_frac(row,2)));
    // tripletList.push_back(Eigen::Triplet<double>(row, col+1, (1-P_frac(row,0)) * P_frac(row,1) * P_frac(row,2)));
    // tripletList.push_back(Eigen::Triplet<double>(row, col+nx, P_frac(row,0) * (1-P_frac(row,1)) * P_frac(row,2)));
    // tripletList.push_back(Eigen::Triplet<double>(row, col+nx*ny, P_frac(row,0) * P_frac(row,1) * (1-P_frac(row,2))));
    // tripletList.push_back(Eigen::Triplet<double>(row, col+1+nx, (1-P_frac(row,0)) * (1-P_frac(row,1)) * P_frac(row,2)));
    // tripletList.push_back(Eigen::Triplet<double>(row, col+1+nx*ny, (1-P_frac(row,0)) * P_frac(row,1) * (1-P_frac(row,2))));
    // tripletList.push_back(Eigen::Triplet<double>(row, col+nx+nx*ny, (P_frac(row,0) * (1-P_frac(row,1)) * (1-P_frac(row,2)))));
    // tripletList.push_back(Eigen::Triplet<double>(row, col+1+nx+nx*ny, (1-P_frac(row,0)) * (1-P_frac(row,1)) * (1-P_frac(row,2))));
    

    tripletList.push_back(Eigen::Triplet<double>(row, col, (1-w_x) * (1-w_y) * (1-w_z)));
    tripletList.push_back(Eigen::Triplet<double>(row, col+1, w_x * (1-w_y) * (1-w_z)));
    tripletList.push_back(Eigen::Triplet<double>(row, col+nx, (1-w_x) * w_y * (1-w_z)));
    tripletList.push_back(Eigen::Triplet<double>(row, col+nx*ny, (1-w_x) * (1-w_y) * w_z));
    tripletList.push_back(Eigen::Triplet<double>(row, col+1+nx, w_x * w_y * (1-w_z)));
    tripletList.push_back(Eigen::Triplet<double>(row, col+1+nx*ny, w_x * (1-w_y) * w_z));
    tripletList.push_back(Eigen::Triplet<double>(row, col+nx+nx*ny, (1-w_x) * w_y * w_z));
    tripletList.push_back(Eigen::Triplet<double>(row, col+1+nx+nx*ny, w_x * w_y * w_z));
    
  }
  W.setFromTriplets(tripletList.begin(), tripletList.end());
  ////////////////////////////////////////////////////////////////////////////
}
