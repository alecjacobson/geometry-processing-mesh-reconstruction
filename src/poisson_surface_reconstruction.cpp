#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <iostream>
#include "fd_interpolate.h"
#include "fd_partial_derivative.h"
#include "fd_grad.h"

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
  // distribute the given normals N 
  Eigen::SparseMatrix<double> Wx, Wy, Wz;

  // Moving P to the left by h/2 is the same as moving corner to the right by h/2
  fd_interpolate((nx - 1), ny, nz, h, corner + Eigen::RowVector3d(-h/2, 0, 0), P, Wx);
  fd_interpolate(nx, (ny - 1), nz, h, corner + Eigen::RowVector3d(0, -h/2, 0), P, Wy);
  fd_interpolate(nx, ny, (nz - 1), h, corner + Eigen::RowVector3d(0, 0, -h/2), P, Wz);

  Eigen::VectorXd Vx(Wx.cols());
  Vx << Wx.transpose() * N.col(0);
  Eigen::VectorXd Vy(Wy.cols());
  Vy << Wy.transpose() * N.col(1);
  Eigen::VectorXd Vz(Wz.cols());
  Vz << Wz.transpose() * N.col(2);
  Eigen::VectorXd my_vec(Vx.size() + Vy.size() + Vz.size());
  my_vec << Vx, Vy, Vz;


  // Compute Gradient:
  Eigen::SparseMatrix<double> G;
  fd_grad(nx, ny, nz, h, G);


  // Compute g: G.T G g = G.T v
  // Idea inspired by BiCGSTAB library
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
  cg.compute(G.transpose() * G);
  g = cg.solve(G.transpose() * my_vec);

  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // Find W on the primary (non-staggered) grid:
  Eigen::SparseMatrix<double> W;
  fd_interpolate(nx, ny, nz, h, corner, P, W);

  // Compute sigma:
  double sigma = (W * g).sum() / W.rows();
  g = g - Eigen::VectorXd::Ones(g.size()) * sigma;

  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
