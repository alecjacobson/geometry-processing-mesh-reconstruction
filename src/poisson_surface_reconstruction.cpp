#include "poisson_surface_reconstruction.h"
#include "fd_grad.h"
#include "fd_interpolate.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>

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

  // Get staggered gradient matrix
  Eigen::SparseMatrix<double> G((nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1), nx*ny*nz);
  fd_grad(nx, ny, nz, h, G);

  // Get interpolated normals
  Eigen::SparseMatrix<double> Wx(n, (nx-1)*ny*nz), Wy(n, nx*(ny-1)*nz), Wz(n, nx*ny*(nz-1));

  // Get Wx
  corner(0) = corner(0) + h/2;
  fd_interpolate(nx-1, ny, nz, h, corner, P, Wx);
  corner(0) = corner(0) - h/2;  

  // Get Wy
  corner(1) = corner(1) + h/2;
  fd_interpolate(nx, ny-1, nz, h, corner, P, Wy);
  corner(1) = corner(1) - h/2;

  // Get Wz
  corner(2) = corner(2) + h/2;
  fd_interpolate(nx, ny, nz-1, h, corner, P, Wz);
  corner(2) = corner(2) - h/2;  

  // Project normals to staggered grids
  Eigen::VectorXd vx, vy, vz, tmp, v(3*nx*ny*nz - nx*ny - ny*nz - nx*nz);
  vx = Eigen::SparseMatrix<double>(Wx.transpose())*N.col(0);
  vy = Eigen::SparseMatrix<double>(Wy.transpose())*N.col(1);
  vz = Eigen::SparseMatrix<double>(Wz.transpose())*N.col(2);

  igl::cat(1, vx, vy, tmp);
  igl::cat(1, tmp, vz, v);

  // Solve
  Eigen::SparseMatrix<double> Gt = Eigen::SparseMatrix<double>(G.transpose());
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(Gt*G);
  g = solver.solve(Gt*v);


  // Get interpolation matrix
  Eigen::SparseMatrix<double> W(n, nx*ny*nz);
  fd_interpolate(nx, ny, nz, h, corner, P, W);
  Eigen::MatrixXd ones = (Eigen::MatrixXd::Zero(1, n).array() + 1).matrix(); 

  // Get sigma
  Eigen::VectorXd sigma = (ones*W*g)/n;
  
  // Shift
  g = (g.array() - sigma(0,0)).matrix();

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}