#ifndef FD_GRAD_H
#define FD_GRAD_H
#include <Eigen/Sparse>
// Construct a gradient matrix for a finite-difference grid
//
// Inputs:
//   nx  number of grid steps along the x-direction
//   ny  number of grid steps along the y-direction
//   nz  number of grid steps along the z-direction
//   h  grid step size
// Outputs:
//   G  (nx-1)*ny*nz+ nx*(ny-1)*nz+ nx*ny*(nz-1) by nx*ny*nz sparse gradient
//     matrix: G = [Dx;Dy;Dz]
//
// See also: fd_partial_derivative.h
void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G);
#endif
