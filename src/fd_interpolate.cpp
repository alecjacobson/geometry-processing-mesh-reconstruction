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
  // set size of W
  W.resize(P.rows(), nx*ny*nz);
  int pi, pj, pk, col;
  
  for(int p = 0; p < P.rows(); p++) {
    
  }
}
