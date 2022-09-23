#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include "fd_grad.h"
#include "fd_interpolate.h"
#include <iostream>

using namespace std;
using namespace Eigen;

void poisson_surface_reconstruction(
    const MatrixXd & P,
    const MatrixXd & N,
    MatrixXd & V,
    MatrixXi & F)
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
  RowVector3d corner = P.colwise().minCoeff().array()-pad*h;
  // Grid dimensions should be at least 3
  nx = max((P.col(0).maxCoeff()-P.col(0).minCoeff()+(2.*pad)*h)/h,3.);
  ny = max((P.col(1).maxCoeff()-P.col(1).minCoeff()+(2.*pad)*h)/h,3.);
  nz = max((P.col(2).maxCoeff()-P.col(2).minCoeff()+(2.*pad)*h)/h,3.);
  // Compute positions of grid nodes
  MatrixXd x(nx*ny*nz, 3);
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < ny; j++)
    {
      for(int k = 0; k < nz; k++)
      {
         // Convert subscript to index
         const auto ind = i + nx*(j + k * ny);
         x.row(ind) = corner + h * RowVector3d(i,j,k);
      }
    }
  }
  VectorXd g = VectorXd::Zero(nx*ny*nz);

  VectorXd vx, vy, vz, v;
  SparseMatrix<double> Wx, Wy, Wz;
  //notice the shift of corner here
  fd_interpolate(nx - 1, ny, nz, h, corner + RowVector3d(0.5 * h, 0, 0), P, Wx);
  vx = Wx.transpose() * N.col(0);
  fd_interpolate(nx, ny - 1, nz, h, corner + RowVector3d(0, 0.5 * h, 0), P, Wy);
  vy = Wy.transpose() * N.col(1);
  fd_interpolate(nx, ny, nz - 1, h, corner + RowVector3d(0, 0, 0.5 * h), P, Wz);
  vz = Wz.transpose() * N.col(2);

  v.resize(vx.rows() + vy.rows() + vz.rows()); // remember to initialize size of v here
  v << vx, vy, vz;

  SparseMatrix<double> G, left, W; // product of sparsematrices are also sparse
  VectorXd right;
  fd_grad(nx, ny, nz, h, G);
  left = G.transpose() * G;
  right = G.transpose() * v;
  BiCGSTAB<SparseMatrix<double>> solver;
  solver.compute(left);
  g = solver.solve(right);

  fd_interpolate(nx, ny, nz, h, corner, P, W);

  // compute sigma and then shift g
  g = (g.array() - (MatrixXd::Ones(1, n) * W * g).value() * 1.0 / n).matrix();


  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
