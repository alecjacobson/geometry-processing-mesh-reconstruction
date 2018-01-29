#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <iostream>
#include <igl/cat.h>

using namespace std;

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  Eigen::SparseMatrix<double> D_x((nx-1) * ny * nz, nx * ny * nz);
  fd_partial_derivative(nx, ny, nz, h, 0, D_x);
  cout << "Done D_x " << D_x.rows() << ' ' << D_x.cols() << endl;

  Eigen::SparseMatrix<double> D_y(nx * (ny-1) * nz, nx * ny * nz);
  fd_partial_derivative(nx, ny, nz, h, 1, D_y);
  cout << "Done D_y " << D_y.rows() << ' ' << D_y.cols() << endl;

  Eigen::SparseMatrix<double> D_z(nx * ny * (nz-1), nx * ny * nz);
  fd_partial_derivative(nx, ny, nz, h, 2, D_z);
  cout << "Done D_z " << D_z.rows() << ' ' << D_z.cols() << endl;

  // Concatenate using igl concatenation
  Eigen::SparseMatrix<double> tempMatrix(D_x.rows() + D_y.rows(), G.cols());
  tempMatrix.reserve(D_x.nonZeros() + D_y.nonZeros());
  igl::cat(1, D_x, D_y, tempMatrix);

  G.reserve(D_x.nonZeros() + D_y.nonZeros() + D_z.nonZeros());
  igl::cat(1, tempMatrix, D_z, G);
}
