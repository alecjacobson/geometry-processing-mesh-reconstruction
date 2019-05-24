#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <igl/cat.h>
#include <algorithm>
#include "fd_interpolate.h"
#include "fd_grad.h"
#include "fd_get_ind.h"
#include <iostream>

using namespace std;

void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
{
  ////////////////////////////////////////////////////////////////////////////
  // Construct FD grid
  ////////////////////////////////////////////////////////////////////////////
  // number of input points
  const int n = P.rows();
  // Grid dimensions
  int nx, ny, nz;
  // Minimum extent (side length of bounding box) of points
  double min_extent = (P.colwise().maxCoeff()-P.colwise().minCoeff()).minCoeff();
  // padding: number of cells beyond bounding box of input points
  const double pad = 8;
  // choose grid spacing (h) so that shortest side gets 40+2*pad samples
  // 30 would cause some distortion between fingers
  double h  = min_extent/double(40-1);
  // Place bottom-left-front corner of grid at minimum of points minus padding
  Eigen::RowVector3d corner = P.colwise().minCoeff().array()-pad*h;
  // Grid dimensions should be at least 3 
  nx = std::max((P.col(0).maxCoeff()-P.col(0).minCoeff()+(2.*pad)*h)/h,3.);
  ny = std::max((P.col(1).maxCoeff()-P.col(1).minCoeff()+(2.*pad)*h)/h,3.);
  nz = std::max((P.col(2).maxCoeff()-P.col(2).minCoeff()+(2.*pad)*h)/h,3.);
  cout << nx << ' ' << ny << ' ' << nz << endl;
  // Compute positions of grid nodes
  Eigen::MatrixXd x(nx*ny*nz, 3);
  for(int i = 0; i < nx; i++) 
  {
    for(int j = 0; j < ny; j++)
    {
      for(int k = 0; k < nz; k++)
      {
         // Convert subscript to index
         const auto ind = fd_get_ind(nx, ny, nz, i, j, k); 
         x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  int grid_volume = nx * ny * nz;
  int x_grid_volume = (nx - 1) * ny * nz;
  int y_grid_volume = nx * (ny - 1) * nz;
  int z_grid_volume = nx * ny * (nz - 1);
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // Perform Reverse Interpolation onto all staggered grid coordinates (in x, y, and z)
  ////////////////////////////////////////////////////////////////////////////
  Eigen::SparseMatrix<double> W(n, grid_volume);
  fd_interpolate(nx, ny, nz, h, corner, P, W);

  Eigen::SparseMatrix<double> W_x(n, x_grid_volume), W_y(n, y_grid_volume), W_z(n, z_grid_volume);
  fd_interpolate(nx, ny, nz, h, corner + h/2.0 * Eigen::RowVector3d(1,0,0), P, W_x);
  fd_interpolate(nx, ny, nz, h, corner + h/2.0 * Eigen::RowVector3d(0,1,0), P, W_y);
  fd_interpolate(nx, ny, nz, h, corner + h/2.0 * Eigen::RowVector3d(0,0,1), P, W_z);

  ////////////////////////////////////////////////////////////////////////////
  // Perform Reverse Interpolation onto all grid coordinates
  ////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXd N_x = N.col(0);
  Eigen::VectorXd v_x = W_x.transpose() * N_x;
  cout << "v_x shape: " << v_x.rows() << ' ' << v_x.cols() << endl;

  Eigen::VectorXd N_y = N.col(1);
  Eigen::VectorXd v_y = W_y.transpose() * N_y;
  cout << "v_y shape: " << v_y.rows() << ' ' << v_y.cols() << endl;

  Eigen::VectorXd N_z = N.col(2);
  Eigen::VectorXd v_z = W_z.transpose() * N_z;
  cout << "v_z shape: " << v_z.rows() << ' ' << v_z.cols() << endl;

  Eigen::VectorXd v(x_grid_volume + y_grid_volume + z_grid_volume);
  v << v_x, v_y, v_z;
  cout << "v_stacked shape: " << v.rows() << ' ' << v.cols() << endl;
  
  ////////////////////////////////////////////////////////////////////////////
  // Get Sparse gradient matrix
  ////////////////////////////////////////////////////////////////////////////
  Eigen::SparseMatrix<double> G(x_grid_volume + y_grid_volume + z_grid_volume, grid_volume);
  fd_grad(nx, ny, nz, h, G);
  cout << "G shape: " << G.rows() << ' ' << G.cols() << endl;
  
  ////////////////////////////////////////////////////////////////////////////
  // Optimize objective
  ////////////////////////////////////////////////////////////////////////////
  Eigen::SparseMatrix<double> G_transpose = G.transpose();
  
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > cg_solver;
  cg_solver.compute(G_transpose * G);
  
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);
  g = cg_solver.solve(G_transpose * v);
  cout << "g shape: " << g.rows() << ' ' << g.cols() << endl;
  cout << "g min max avg: " << g.minCoeff() << ' ' << g.maxCoeff() << ' ' << g.sum() / (1.0*n) << endl;
  
  ////////////////////////////////////////////////////////////////////////////
  // Interpolate to find appropriate sigma surface
  ////////////////////////////////////////////////////////////////////////////
  Eigen::VectorXd sigma_calculation = Eigen::MatrixXd::Ones(1, n) * W * g;
  cout << "sigma shape: " << sigma_calculation.rows() << ' ' << sigma_calculation.cols() << endl;
  double sigma = sigma_calculation.value() / n;
  cout << "sigma value: " << sigma << endl;

  // Subtract such that when g = sigma, the result is zero
  g = g.array() - sigma;
  cout << "adjusted g SHAPE: " << g.rows() << ' ' << g.cols() << endl;    

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
