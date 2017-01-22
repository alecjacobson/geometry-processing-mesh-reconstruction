#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <fd_interpolate.h>
#include <fd_grad.h>
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
  Eigen::SparseMatrix<double> Wx(n,(nx-1)*ny*nz);
  Eigen::RowVector3d cornerDx = corner;
  cornerDx(0) = cornerDx(0) + h;
  Eigen::MatrixXd Nx = N;
  Nx.col(1) = Eigen::MatrixXd::Zero(Nx.rows(), 1);
  Nx.col(2) = Eigen::MatrixXd::Zero(Nx.rows(), 1);
  fd_interpolate(nx, ny, nz, h, cornerDx, Nx, Wx);
  Eigen::MatrixXd vx = Wx.transpose()*Nx;

  Eigen::SparseMatrix<double> Wy(n,nx*(ny-1)*nz);
  Eigen::RowVector3d cornerDy = corner;
  cornerDy(1) = cornerDy(1) + h;
  Eigen::MatrixXd Ny = N;
  Ny.col(0) = Eigen::MatrixXd::Zero(Ny.rows(), 1);
  Ny.col(2) = Eigen::MatrixXd::Zero(Ny.rows(), 1);
  fd_interpolate(nx, ny, nz, h, cornerDy, Ny, Wy);
  Eigen::MatrixXd vy = Wy.transpose()*Ny;

  Eigen::SparseMatrix<double> Wz(n, nx*ny*(nz-1));
  Eigen::RowVector3d cornerDz = corner;
  cornerDz(2) = cornerDz(2) + h;
  Eigen::MatrixXd Nz = N;
  Nz.col(0) = Eigen::MatrixXd::Zero(Nz.rows(), 1);
  Nz.col(1) = Eigen::MatrixXd::Zero(Nz.rows(), 1);
  fd_interpolate(nx, ny, nz, h, cornerDz, Nz, Wz);
  Eigen::MatrixXd vz = Wz.transpose()*Nz;

  Eigen::MatrixXd v = Eigen::MatrixXd(vx.rows() + vy.rows() + vz.rows(), 1);
  v.block(0, 0, vx.rows(), 1) = vx.col(0);
  v.block(vx.rows(), 0, vy.rows(), 1) = vy.col(1);
  v.block(vx.rows() + vy.rows(), 0, vz.rows(), 1) = vz.col(2);

  Eigen::SparseMatrix<double> G((nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1),nx*ny*nz);
  fd_grad(nx, ny, nz, h, G);

  Eigen::SparseMatrix<double> A = G.transpose()*G;
  Eigen::MatrixXd b = G.transpose()*v;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  g = solver.solve(b);

  
  //Eigen::SparseMatrix<double> WxT = Wx.transpose();
  //Eigen::SparseMatrix<double> WyT = Wy.transpose();
  //Eigen::SparseMatrix<double> WzT = Wz.transpose();
  //Eigen::SparseMatrix<double> W(WxT.rows() + WyT.rows() + WzT.rows(),n);
  //int nnz = WxT.nonZeros() + WyT.nonZeros() + WzT.nonZeros();
  //W.reserve(nnz);
  //std::vector<Eigen::Triplet<double> > tripletList;
  //tripletList.reserve(nnz);
  //for (int c = 0; c < n; ++c) {
	 // for (Eigen::SparseMatrix<double>::InnerIterator it(WxT, c); it; ++it) {
		//  tripletList.push_back(Eigen::Triplet<double>(it.row(), c, it.value()));
	 // }
	 // for (Eigen::SparseMatrix<double>::InnerIterator it(WyT, c); it; ++it) {
		//  tripletList.push_back(Eigen::Triplet<double>(it.row(), c, it.value()));
	 // }
	 // for (Eigen::SparseMatrix<double>::InnerIterator it(WzT, c); it; ++it) {
		//  tripletList.push_back(Eigen::Triplet<double>(it.row(), c, it.value()));
	 // }
  //}
  //W.setFromTriplets(tripletList.begin(), tripletList.end());
  //W.makeCompressed();
  //W = W.transpose();

  Eigen::SparseMatrix<double> Wg(n, nx*ny*nz);
  fd_interpolate(nx, ny, nz, h, corner, g, Wg);

  Eigen::MatrixXd sigma = (1 / n)*(Eigen::MatrixXd::Ones(1, n))*Wg*g;

  g = g - sigma;

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
