#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include "fd_interpolate.h"
#include "fd_grad.h"
#include <algorithm>
#include <iostream>
#include <Eigen/Sparse>


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

  /* Probably, the definition of h should be instead as follows to ensure atleast 30+2*pad samples along the *smallest* dimension:
        double min_extent = (P.colwise().maxCoeff()-P.colwise().minCoeff()).minCoeff();
        double h = min_extent/double(30+2*pad);
  */

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

  // to hold the weight matrix for staggered interpolation
  Eigen::SparseMatrix<double> Wx;
  Eigen::SparseMatrix<double> Wy;
  Eigen::SparseMatrix<double> Wz; 


  // need to interpolate only the staggered positions
  fd_interpolate(nx-1, ny, nz, h, corner + (h/2)*Eigen::RowVector3d(1,0,0), P, Wx);
  fd_interpolate(nx, ny-1, nz, h, corner + (h/2)*Eigen::RowVector3d(0,1,0), P, Wy);
  fd_interpolate(nx, ny, nz-1, h, corner + (h/2)*Eigen::RowVector3d(0,0,1), P, Wz);

  
  Eigen::MatrixXd N_interpolated_x = Wx.transpose() * N.col(0);
  Eigen::MatrixXd N_interpolated_y = Wy.transpose() * N.col(1);
  Eigen::MatrixXd N_interpolated_z = Wz.transpose() * N.col(2);

  // concatenate
  Eigen::MatrixXd N_interpolated(N_interpolated_x.rows()+N_interpolated_y.rows()+N_interpolated_z.rows(), N_interpolated_x.cols());

  N_interpolated << N_interpolated_x, N_interpolated_y, N_interpolated_z;




  // compute the gradient
  Eigen::SparseMatrix<double> G;
  
  fd_grad(nx, ny, nz, h, G);
  

  // solve the system
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  solver.compute(G.transpose()*G);

  Eigen::VectorXd g_without_sigma = Eigen::VectorXd::Zero(nx*ny*nz);

  g_without_sigma = solver.solve(G.transpose()*N_interpolated);

  // Shift by the iso-value
  Eigen::VectorXd n_ones = Eigen::VectorXd::Ones(n);


  // to hold the weight matrix for regular grid (g) interpolation
  Eigen::SparseMatrix<double> W;
  fd_interpolate(nx, ny, nz, h, corner, P, W);

  double sigma;
  sigma = (1./n)*n_ones.transpose()*W*g;

  g = g_without_sigma - sigma*Eigen::VectorXd::Ones(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}

