#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <fd_interpolate.h>
#include <fd_grad.h>
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
  Eigen::MatrixXd grid(nx*ny*nz, 3);
  for(int i = 0; i < nx; i++) 
  {
    for(int j = 0; j < ny; j++)
    {
      for(int k = 0; k < nz; k++)
      {
         // Convert subscript to index
         const auto ind = i + nx*(j + k * ny);
         grid.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
      }
    }
  }
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  // Compute vx,vy,vz on staggered grid
  Eigen::SparseMatrix<double> Wx(n,(nx-1)*ny*nz);
  fd_interpolate(nx-1, ny, nz, h, corner + 0.5*h*Eigen::RowVector3d(1,0,0), P, Wx);
  Eigen::MatrixXd vx = Wx.transpose()*N.col(0);

  Eigen::SparseMatrix<double> Wy(n,nx*(ny-1)*nz);
  fd_interpolate(nx, ny-1, nz, h, corner + 0.5*h*Eigen::RowVector3d(0, 1, 0), P, Wy);
  Eigen::MatrixXd vy = Wy.transpose()*N.col(1);

  Eigen::SparseMatrix<double> Wz(n, nx*ny*(nz-1));
  fd_interpolate(nx, ny, nz-1, h, corner + 0.5*h*Eigen::RowVector3d(0, 0, 1), P, Wz);
  Eigen::MatrixXd vz = Wz.transpose()*N.col(2);

  // Combine vx,vy,vz into v
  Eigen::MatrixXd v = Eigen::MatrixXd(vx.rows() + vy.rows() + vz.rows(), 1);
  v.block(0, 0, vx.rows(), 1) = vx;
  v.block(vx.rows(), 0, vy.rows(), 1) = vy;
  v.block(vx.rows() + vy.rows(), 0, vz.rows(), 1) = vz;

  // Compute grad matrix G
  Eigen::SparseMatrix<double> G((nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1),nx*ny*nz);
  fd_grad(nx, ny, nz, h, G);

  // Solve linear system for g
  Eigen::SparseMatrix<double> A = G.transpose()*G;
  Eigen::MatrixXd b = G.transpose()*v;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  g = solver.solve(b);

  // Compute weight matrix W for primary grid
  Eigen::SparseMatrix<double> W(n,nx*ny*nz);
  fd_interpolate(nx, ny, nz, h, corner, P, W);

  // Average results for points in cloud
  double sigma = (W*g).sum()/n;

  // Set up isosurface
  g = g - sigma*Eigen::MatrixXd::Ones(nx*ny*nz,1);

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, grid, nx, ny, nz, V, F);

}
