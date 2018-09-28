#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include "fd_grad.h"
#include "fd_interpolate.h"

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
  ////////////////////////////////////////////////////////////////////////////
  //Interpolation matrices on the staggered grids and the primary grid
  Eigen::SparseMatrix<double> Wx(P.rows(), (nx - 1) * ny * nz);
  Eigen::SparseMatrix<double> Wy(P.rows(), (ny - 1) * nx * nz);
  Eigen::SparseMatrix<double> Wz(P.rows(), (nz - 1) * nx * ny);
  Eigen::SparseMatrix<double> W(P.rows(), nx*ny*nz);
  Eigen::SparseMatrix<double> G;
  // Corner points of the staggered grids
  Eigen::RowVector3d cornerx, cornery, cornerz;
  cornerx(0) = corner(0) + h/2;  cornerx(1) = corner(1);  cornerx(2) = corner(2);
  cornery(0) = corner(0);  cornery(1) = corner(1) + h/2;  cornery(2) = corner(2);
  cornerz(0) = corner(0);  cornerz(1) = corner(1);  cornerz(2) = corner(2) + h/2;

  // Calculate interpolation matrix on the staggered grids
  fd_interpolate(nx - 1, ny, nz, h, cornerx, P, Wx);
  fd_interpolate(nx, ny - 1, nz, h, cornery, P, Wy);
  fd_interpolate(nx, ny, nz - 1, h, cornerz, P, Wz);

  // Calculate interpolation matrix W on the primary grid
  fd_interpolate(nx, ny, nz, h, corner, P, W);
  Eigen::MatrixXd vx(P.rows(), 1), vy(P.rows(), 1), vz(P.rows(), 1);
  vx = Wx.transpose() * N.col(0);
  vy = Wy.transpose() * N.col(1);
  vz = Wz.transpose() * N.col(2);

  // Concatenate vx, vy and vz
  Eigen::VectorXd v(vx.rows()+vy.rows()+vz.rows());
  for (int i = 0; i < vx.rows(); i++) {
	  v(i) = vx(i);
  }
  for (int i = 0; i < vy.rows(); i++) {
	  v(i+vx.rows()) = vy(i);
  }
  for (int i = 0; i < vz.rows(); i++) {
	  v(i+vx.rows()+vy.rows()) = vz(i);
  }
 
  fd_grad(nx, ny, nz, h, G);

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
  g = solver.compute(G.transpose() * G).solve(G.transpose() * v);

  double sigma = (1/(double) P.rows()) * (W * g).sum();

  Eigen::VectorXd ones(g.size());
  ones.setOnes();
  g = g - sigma * ones;

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
