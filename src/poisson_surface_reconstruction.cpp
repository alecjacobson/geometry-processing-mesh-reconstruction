#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <igl/cat.h>
#include "fd_interpolate.h"
#include "fd_grad.h"
#include <fstream>

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
  //printf("%d %d %d %d ",n,nx,ny,nz);
  Eigen::RowVector3d c1(h/2,0,0),c2(0,h/2,0),c3(0,0,h/2);
  Eigen::SparseMatrix<double> W1,W2,W3,W;
  fd_interpolate(nx-1,ny,nz,h,corner+c1,P,W1);
  fd_interpolate(nx,ny-1,nz,h,corner+c2,P,W2);
  fd_interpolate(nx,ny,nz-1,h,corner+c3,P,W3);
  Eigen::VectorXd V1((nx-1)*ny*nz),V2(nx*(ny-1)*nz),V3(nx*ny*(nz-1)),vv((nx-1)*ny*nz+nx*(ny-1)*nz+nx*ny*(nz-1));
  V1=W1.transpose()*N.col(0);
  V2=W2.transpose()*N.col(1);
  V3=W3.transpose()*N.col(2);
  vv << V1,V2,V3;

  Eigen::SparseMatrix<double> gg;
  fd_grad(nx,ny,nz,h,gg);

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(gg.transpose()*gg);
  printf("%d %d ",gg.rows(),gg.cols());
  Eigen::MatrixXd b=gg.transpose()*vv;
  g=solver.solve(b);

  fd_interpolate(nx,ny,nz,h,corner,P,W);
  double sigma=(Eigen::MatrixXd::Ones(1,n)*W*g).value()/n;
  g=g-Eigen::VectorXd::Ones(nx*ny*nz)*sigma;
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
