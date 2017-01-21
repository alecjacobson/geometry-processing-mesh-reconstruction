#ifndef FD_PARTIAL_DERIVATIVE_H
#define FD_PARTIAL_DERIVATIVE_H
#include <Eigen/Sparse>
// Construct a partial derivative matrix for a finite-difference grid in a
// given direction. Derivative are computed using first-order differences onto
// a staggered grid
//
// Inputs:
//   nx  number of grid steps along the x-direction
//   ny  number of grid steps along the y-direction
//   nz  number of grid steps along the z-direction
//   h  grid step size
//   dir  index indicating direction: 0-->x, 1-->y, 2-->z
// Outputs:
//   D  m by nx*ny*nz sparse partial derivative matrix, where:
//     m = (nx-1)*ny*nz  if dir = 0
//     m = nx*(ny-1)*nz  if dir = 1
//     m = nx*ny*(nz-1)  otherwise (if dir = 2)
//
// See also: fd_partial_derivative.h
void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D);
#endif
