#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <iostream>
#include <Eigen/Sparse>
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
  ////////////////////////////////////////////////////////////////////////////
  Eigen::SparseMatrix<double> Dx((nx - 1) * ny * nz, nx * ny * nz);
  fd_partial_derivative(nx, ny, nz, h, 0, Dx);  
  std::cout << "Dx Done" << std::endl;  
  Eigen::SparseMatrix<double> Dy(nx * (ny - 1) * nz, nx * ny * nz);
  fd_partial_derivative(nx, ny, nz, h, 1, Dy);  
  std::cout << "Dy Done" << std::endl;
  Eigen::SparseMatrix<double> Dz(nx * ny * (nz-1), nx * ny * nz);
  fd_partial_derivative(nx, ny, nz, h, 2, Dz);  
  std::cout << "Dz Done" << std::endl;

  Eigen::SparseMatrix<double> G_temp((nx - 1) * ny * nz + nx * (ny -1) * nz, nx * ny * nz);
  igl::cat(1, Dx, Dy, G_temp);
  igl::cat(1, G_temp, Dz, G); 
}



