#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
/////////////////////////////////////////
//My Code
#include <Eigen/Sparse>
#include <iostream>
#include "fd_interpolate.h"
#include "fd_grad.h"

//End of My Code
//////////////////////////////////////////
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
    for(int j = 0; j < ny; j++)
      for(int k = 0; k < nz; k++)
      {
         // Convert subscript to index
         const auto ind = i + nx*(j + k * ny);
         x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
      }
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  // Compute grad matrix G
  Eigen::SparseMatrix<double> G;
  fd_grad(nx, ny, nz, h, G);

  std::vector<Eigen::SparseMatrix<double>> W_xyz(3);
  std::vector<Eigen::VectorXd> v_xyz(3);
  for (int i = 0; i < W_xyz.size(); i++)
  {
    double offset[3] = { 0 };
    offset[i] = 1;
    W_xyz.at(i).resize(n, (nx - offset[0])*(ny - offset[1])*(nz - offset[2]));
    fd_interpolate((nx - offset[0]), (ny - offset[1]), (nz - offset[2]), h, corner, P, W_xyz.at(i));
    v_xyz.at(i)= W_xyz.at(i).transpose()*N.col(i);
  }

  //from https://stackoverflow.com/questions/25691324/how-to-concatenate-vectors-in-eigen
  Eigen::VectorXd v(v_xyz[0].rows()+ v_xyz[1].rows()+ v_xyz[2].rows());
  v << v_xyz[0], v_xyz[1], v_xyz[2];

  Eigen::SparseMatrix<double> GG = G.transpose()*G;
  Eigen::MatrixXd Gv = G.transpose()*v;

  //from https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient.html
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::UpLoType::Lower | Eigen::UpLoType::Upper> cg;
  cg.compute(GG);
  g = cg.solve(Gv);
  std::cout << "#iterations:     " << cg.iterations() << std::endl;
  std::cout << "estimated error: " << cg.error() << std::endl;

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}

