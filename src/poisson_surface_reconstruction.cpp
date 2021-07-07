#include "poisson_surface_reconstruction.h"
#include "fd_interpolate.h"
#include "fd_grad.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <iostream>

using namespace std;
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
  
  // debug W
  // Eigen::SparseMatrix<double> W;
  // fd_interpolate(nx, ny, nz, h, corner, P, W);
  // for (int ii = 0; ii < 100; ii++){
  //   cout << P.row(ii) << ",";
  //   cout << (W*x).row(ii) << endl;
  // }

  // construct Wx Wy Wz
  // cout << "start interpolation ...";
  Eigen::SparseMatrix<double> Wx(n, (nx-1)*ny*nz), Wy(n, nx*(ny-1)*nz), Wz(n, nx*ny*(nz-1)); 
  Eigen::RowVector3d tempVec(3);
  tempVec << 0.5*h, 0, 0;
  fd_interpolate(nx-1, ny, nz, h, corner+tempVec, P, Wx);
  tempVec << 0, 0.5*h, 0;
  fd_interpolate(nx, ny-1, nz, h, corner+tempVec, P, Wy);
  tempVec << 0, 0, 0.5*h;
  fd_interpolate(nx, ny, nz-1, h, corner+tempVec, P, Wz);
  // cout << "finish" << endl;

  // construct v
  // cout << "construct v ...";
  Eigen::VectorXd v(Wx.cols() + Wy.cols() + Wz.cols());
  Eigen::VectorXd vx = Wx.transpose() * N.col(0);
  Eigen::VectorXd vy = Wy.transpose() * N.col(1);
  Eigen::VectorXd vz = Wz.transpose() * N.col(2);

  // cout << vx.cols() << "," << vy.cols() << "," << vz.cols() << endl;
  for (int ii = 0; ii < Wx.cols(); ii++){
    v(ii) = vx(ii);
  }
  for (int ii = 0; ii < Wy.cols(); ii++){
    v(Wx.cols() + ii) = vy(ii);
  }
  for (int ii = 0; ii < Wz.cols(); ii++){
    v(Wx.cols() + Wy.cols() + ii) = vz(ii);
  }
  // cout << "finish" << endl;

  // construct G
  // cout << "construct G ...";
  Eigen::SparseMatrix<double>  G((nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1), nx*ny*nz);
  fd_grad(nx, ny, nz, h, G);
  // cout << "finish" << endl;

  // solve 
  // cout << "solve ...";
  Eigen::SparseMatrix<double> A = G.transpose() * G;
  Eigen::VectorXd b = G.transpose() * v;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > conjGradSolve;
  conjGradSolve.compute(A);
  g = conjGradSolve.solve(b);
  // cout << "finish" << endl;

  // determine sigma
  Eigen::SparseMatrix<double> W;
  fd_interpolate(nx, ny, nz, h, corner, P, W);
  auto mat = (W * g);
  double sigma = mat.array().sum() / W.rows();
  g.array() -= sigma;

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
