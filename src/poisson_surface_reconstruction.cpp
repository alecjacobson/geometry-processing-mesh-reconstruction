#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include "fd_partial_derivative.h"
#include "fd_grad.h"
#include "fd_interpolate.h" 
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
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
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////

  Eigen::SparseMatrix<double> Weightx(n,(nx-1)*ny*nz);
  Eigen::RowVector3d cornerx = corner;
  cornerx[0] = corner[0] + h/2;
  fd_interpolate(nx-1,ny,nz,h,cornerx,P,Weightx);

  Eigen::SparseMatrix<double> Weighty(n,nx*(ny-1)*nz);
  Eigen::RowVector3d cornery = corner;
  cornery[1] = corner[1] + h/2;
  fd_interpolate(nx,ny-1,nz,h,corner,P,Weighty);

  Eigen::SparseMatrix<double> Weightz(n,nx*ny*(nz-1));
  Eigen::RowVector3d cornerz = corner;
  cornerz[2] = corner[2] + h/2;
  fd_interpolate(nx,ny,nz-1,h,corner,P,Weightz);

  Eigen::SparseMatrix<double> W(n,nx*ny*nz);
  fd_interpolate(nx,ny,nz,h,corner,P, W);

  Eigen::SparseMatrix<double> Gradient;
  fd_grad(nx,ny,nz,h,Gradient);
  
  Eigen::VectorXd vx(n,1), vy(n,1), vz(n,1);
  vx=Weightx.transpose()*N.col(0);
  vy=Weighty.transpose()*N.col(1);
  vz=Weightz.transpose()*N.col(2);
  Eigen::VectorXd B((nx-1)*ny*nz+ nx*(ny-1)*nz+ nx*ny*(nz-1));
  B << vx,vy,vz;
  
  Eigen::SparseMatrix<double> A = Gradient.transpose()*Gradient;
  Eigen::VectorXd b = Gradient.transpose()*B;
  //Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> bicg;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::UpLoType::Lower | Eigen::UpLoType::Upper> bicg;
  g = bicg.compute(A).solve(b);
  
  Eigen::VectorXd i(g.size());
  i.setOnes();
  double sigma = (1/(n*1.0)) * (W*g).sum();
  g = g - sigma*i;
  
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
