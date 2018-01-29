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
  int i_int, j_int, k_int;
  double i, j, k;
  double w_x, w_y, w_z;
  for (int row=0; row<P.rows();i++) {
    i = (P(row,0) - corner(0)) / h;
    j = (P(row,1) - corner(1)) / h;
    k = (P(row,2) - corner(2)) / h;

    i_int = (int)i;
    j_int = (int)j;
    k_int = (int)k;

    w_x = i - (int)i;
    w_y = j - (int)j;
    w_z = k - (int)k;
    const auto col = i_int + nx*(j_int + k_int * ny);
    W.insert(row, col) = (w_x*w_y*w_z);
    W.insert(row, col+1) = ((1-w_x)*w_y*w_z);
    W.insert(row, col+nx) = (w_x*(1-w_y)*w_z);
    W.insert(row, col+nx*ny) = (w_x*w_y*(1-w_z));
    W.insert(row, col+1+nx) = ((1-w_x)*(1-w_y)*w_z);
    W.insert(row, col+1+nx*ny) = ((1-w_x)*w_y*(1-w_z));
    W.insert(row, col+nx*nx+ny) = (w_x*(1-w_y)*(1-w_z));
    W.insert(row, col+1+nx+nx*ny) = ((1-w_x)*(1-w_y)*(1-w_z));
  }
  ////////////////////////////////////////////////////////////////////////////
}
