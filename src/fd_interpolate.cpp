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
  
  // Step over each point in P and find corresponding subscripts
  // on the grid
  for(int p = 0; p < P.rows(); p++) {
    pi = (P(p, 0) - corner(0,0)) / h; // i
    pj = (P(p, 1) - corner(0,1)) / h; // j
    pk = (P(p, 2) - corner(0,2)) / h; // k
    
    // set value in W to indicate node is in P
    col = pi + nx*(pj + pk * ny);
    W.insert(p, col) = 1;
  }
}
