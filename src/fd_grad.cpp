#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <iostream>
#include <igl/cat.h>

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
    //std::cout<<"Grad Here!"<<std::endl;
    Eigen::SparseMatrix<double> Dx( (nx-1) * ny * nz, nx * ny * nz);
    Eigen::SparseMatrix<double> Dy( nx * (ny-1) * nz, nx * ny * nz);
    Eigen::SparseMatrix<double> Dz( nx * ny * (nz-1), nx * ny * nz);
    fd_partial_derivative(nx, ny, nz, h, 0, Dx);
    fd_partial_derivative(nx, ny, nz, h, 1, Dy);
    fd_partial_derivative(nx, ny, nz, h, 2, Dz);

    Eigen::SparseMatrix<double> M(Dx.rows() + Dy.rows(), Dx.cols());
    igl::cat(1, Dx, Dy, M);
    igl::cat(1, M, Dz, G);
    //std::cout<<"Grad Done!"<<std::endl;
  ////////////////////////////////////////////////////////////////////////////
}

