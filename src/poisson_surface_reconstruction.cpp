#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include "fd_grad.h"
#include "fd_interpolate.h"
#include <Eigen/SparseCholesky>
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

  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  Eigen::SparseMatrix<double> W, Wx, Wy, Wz, G;
  fd_interpolate(nx, ny, nz, h, corner, P, W);
  fd_grad(nx, ny, nz, h, G);

  fd_interpolate(nx - 1, ny, nz, h, corner + Eigen::RowVector3d(h * 0.5, 0, 0), P, Wx);
  fd_interpolate(nx, ny - 1, nz, h, corner + Eigen::RowVector3d(0, h * 0.5, 0), P, Wy);
  fd_interpolate(nx, ny, nz - 1, h, corner + Eigen::RowVector3d(0, 0, h * 0.5), P, Wz);

  Eigen::VectorXd vector, vx, vy, vz;
  vx = Wx.transpose() * N.col(0);
  vy = Wy.transpose() * N.col(1);
  vz = Wz.transpose() * N.col(2);

  vector.resize((nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1));
  vector << vx, vy, vz;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  Eigen::SparseMatrix<double> A = G.transpose() * G;
  Eigen::VectorXd b = G.transpose() * vector;
  solver.compute(A);
  if(solver.info()!=Eigen::Success) {
    std::cout << "Decomposition failed" << std::endl;
    return;
  }
  g = solver.solve(b);
  if(solver.info()!=Eigen::Success) {
    std::cout << "Solver failed" << std::endl;
    return;
  }

  double sigma = (1.0 / n) * Eigen::RowVectorXd::Ones(n) * W * g;
  g -= sigma * Eigen::VectorXd::Ones(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
