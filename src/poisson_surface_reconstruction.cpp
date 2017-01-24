#include "poisson_surface_reconstruction.h"
#include "fd_interpolate.h"
#include "fd_grad.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <iostream>

void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
{
  ////////////////////////////////////////////////////////////////////////////
  // Construct FD grid, CONGRATULATIONS! You get this for free!
  ////////////////////////////////////////////////////////////////////////////
  // number of input points
  const int n = P.rows();
  // Grid dimensions
  int nx, ny, nz;
  // Maximum extent (side length of bounding box) of points
  double max_extent =
    (P.colwise().maxCoeff()-P.colwise().minCoeff()).maxCoeff();
  // padding: number of cells beyond bounding box of input points
  const double pad = 8;
  // choose grid spacing (h) so that shortest side gets 30+2*pad samples
  double h  = max_extent/double(30+2*pad);
  // Place bottom-left-front corner of grid at minimum of points minus padding
  Eigen::RowVector3d corner = P.colwise().minCoeff().array()-pad*h;
  // Grid dimensions should be at least 3 
  nx = std::max((P.col(0).maxCoeff()-P.col(0).minCoeff()+(2.*pad)*h)/h,3.);
  ny = std::max((P.col(1).maxCoeff()-P.col(1).minCoeff()+(2.*pad)*h)/h,3.);
  nz = std::max((P.col(2).maxCoeff()-P.col(2).minCoeff()+(2.*pad)*h)/h,3.);
  // Compute positions of grid nodes
  Eigen::MatrixXd x(nx*ny*nz, 3);
  for(int i = 0; i < nx; i++) 
  {
    for(int j = 0; j < ny; j++)
    {
      for(int k = 0; k < nz; k++)
      {
         // Convert subscript to index
         const auto ind = i + nx*(j + k * ny);
         x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
      }
    }
  }
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  int sizeX = (nx-1)*ny*nz;
  int sizeY = nx*(ny-1)*nz;
  int sizeZ = nx*ny*(nz-1);

  // Construct Wx, Wy, Wz
  Eigen::RowVector3d cornerX, cornerY, cornerZ;
  cornerX << corner(0,0) + h / 2, corner(0,1), corner(0,2);
  cornerY << corner(0,0), corner(0,1) + h / 2, corner(0,2);
  cornerZ << corner(0,0), corner(0,1), corner(0,2) + h / 2;
  
  Eigen::SparseMatrix<double> Wx(n, sizeX), Wy(n, sizeY), Wz(n, sizeZ);
  fd_interpolate(nx - 1, ny, nz, h, cornerX, P, Wx);
  fd_interpolate(nx, ny - 1, nz, h, cornerY, P, Wy);
  fd_interpolate(nx, ny, nz - 1, h, cornerZ, P, Wz);
  
  // Construct v
  Eigen::MatrixXd vx(sizeX, 1), vy(sizeY, 1), vz(sizeZ, 1), v(sizeX + sizeY + sizeZ, 1);
  vx = Wx.transpose() * N.col(0);
  vy = Wy.transpose() * N.col(1);
  vz = Wz.transpose() * N.col(2);
  v << vx, vy, vz;
  
  // Construct G
  Eigen::SparseMatrix<double> G(sizeX + sizeY + sizeZ, nx*ny*nz);
  fd_grad(nx, ny, nz, h, G);
  
  // Solve equation
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
  cg.compute(G.transpose() * G);
  g = cg.solve(G.transpose() * v);
  
  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  Eigen::SparseMatrix<double> W(n, nx*ny*nz);
  fd_interpolate(nx, ny, nz, h, corner, P, W);
  double sigma = (1/n) * Eigen::RowVectorXd::Ones(n) * W * g;
  
  g = g - (sigma * Eigen::VectorXd::Ones(nx*ny*nz));
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
