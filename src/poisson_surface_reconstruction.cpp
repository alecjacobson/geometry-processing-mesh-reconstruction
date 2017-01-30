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
  // Construct FD grid, CONGRATULATIONS! You get this for free! -> THANKS BRO!
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
  ////////////////////////////////////////////////////////////////////////////

  fprintf(stderr, "Size of Points: %d x %d\n", P.rows(), P.cols());
  fprintf(stderr, "Grid size: %dx%dx%d with height = %f\n", nx, ny, nz, h);
  fprintf(stderr, "Corder Point: %f, %f, %f\n", corner(0), corner(1), corner(2));  

  //Compute Gradient
  Eigen::SparseMatrix<double> G;
  fd_grad(nx, ny, nz, h, G);  

  //Get interpolation weights for primary and staggered grids (x,y, and z)
  Eigen::SparseMatrix<double> W, Wx, Wy, Wz;  
  fd_interpolate(nx, ny, nz, h, corner, P, W);  	
  fd_interpolate(nx - 1, ny, nz, h, corner + Eigen::RowVector3d(h/2, 0, 0), P, Wx);
  fd_interpolate(nx, ny - 1, nz, h, corner + Eigen::RowVector3d(0, h/2, 0), P, Wy);
  fd_interpolate(nx, ny, nz - 1, h, corner + Eigen::RowVector3d(0, 0, h/2), P, Wz);
  
  //Distrubute the values of the normal to the staggered grid locations, based off the
  //trilinear interpolation weights of the points.	
  Eigen::MatrixXd vx = Wx.transpose() * N.col(0);
  Eigen::MatrixXd vy = Wy.transpose() * N.col(1);
  Eigen::MatrixXd vz = Wz.transpose() * N.col(2);
  Eigen::VectorXd v(Wx.cols() + Wy.cols() + Wz.cols());
  v << vx, vy, vz;  
  
  //Solve the linear system (example from https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient.html)
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
  g = cg.compute(G.transpose() * G).solve(G.transpose() * v);
  std::cout << "#iterations:     " << cg.iterations() << std::endl;
  std::cout << "estimated error: " << cg.error() << std::endl;

  //Interpolate solution g and average at each point
  double sigma = (W * g).array().sum() / P.rows();  
  g.array() -= sigma;


  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
