#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <Eigen/IterativeLinearSolvers>
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

  Eigen::SparseMatrix<double> G, W;
  fd_grad(nx, ny, nz, h, G);
  fd_interpolate(nx, ny, nz, h, corner, P, W);

  Eigen::VectorXd Nvec(N.size());
  Nvec.block(0, 0, N.rows(), 1) = N.col(0);
  Nvec.block(N.rows(), 0, N.rows(), 1) = N.col(1);
  Nvec.block(N.rows() * 2, 0, N.rows(), 1) = N.col(2);
  
  Eigen::VectorXd v = W.transpose()*Nvec;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  Eigen::SparseMatrix<double> A = G.transpose()*G;
  solver.compute(A);

  Eigen::VectorXd b = G.transpose()*v;
  g = solver.solve(b);

  Eigen::RowVectorXd ones = Eigen::RowVectorXd::Ones(P.rows());

  double sigma = ones*W*g;
  g = g - Eigen::VectorXd::Ones(g.rows())*sigma;

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
