#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include "fd_grad.h"
#include "fd_interpolate.h"
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
  //std::cout << "Start" << std::endl;

  //Distribute normal 
  Eigen::VectorXd v = Eigen::VectorXd::Zero((nx-1) * ny * nz + nx * (ny-1) * nz + nx * ny * (nz -1));
  //Vx
  Eigen::SparseMatrix<double> W_x(n, (nx-1) * ny * nz);
  Eigen::VectorXd v_x((nx-1) * ny * nz);
  fd_interpolate(nx-1,ny,nz,h,corner + h * 0.5 * Eigen::RowVector3d(1,0,0), P,W_x);
  v_x = W_x.transpose() * N.col(0);
  //Vy
  Eigen::SparseMatrix<double> W_y(n, nx * (ny-1) * nz);
  Eigen::VectorXd v_y(nx * (ny-1) * nz);
  fd_interpolate(nx,ny-1,nz,h,corner + h * 0.5 * Eigen::RowVector3d(0,1,0), P,W_y);
  v_y = W_y.transpose() * N.col(1);
  //Vz
  Eigen::SparseMatrix<double> W_z(n, nx * ny * (nz-1));
  Eigen::VectorXd v_z(nx * ny * (nz-1));
  fd_interpolate(nx,ny,nz-1,h,corner + h * 0.5 * Eigen::RowVector3d(0,0,1), P,W_z);
  v_z = W_z.transpose() * N.col(2);
  v << v_x, v_y, v_z;
  //std::cout << "Finish normal" << std::endl;
  //solve
  Eigen::SparseMatrix<double> G((nx-1) * ny * nz + nx * (ny-1) * nz + nx * ny * (nz -1), 
    nx * ny * nz);
  fd_grad(nx, ny, nz, h, G);
  //std::cout << "Construct grad" << std::endl;
  
  //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> cg;
  cg.compute(G.transpose() * G);
  Eigen::VectorXd b(nx * ny * nz);
  //std::cout << "Solving" << std::endl;
  b = G.transpose() * v;
  g = cg.solve(b);
  //std::cout << "Solve" << std::endl;
  //sigma
  Eigen::SparseMatrix<double> W_sigma(n, nx * ny * nz);
  fd_interpolate(nx, ny, nz, h, corner, P, W_sigma);
  double sigma = 1.00/n * (Eigen::MatrixXd::Ones(1, n) * W_sigma *g).value();
  //std::cout << "Sigma" << std::endl;
  g = g - sigma * Eigen::MatrixXd::Ones(nx*ny*nz,1);
  //std::cout << "Sigma" << std::endl;
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////

  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
