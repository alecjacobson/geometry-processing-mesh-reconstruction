#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <Eigen/Sparse>
#include "fd_interpolate.h"
#include "fd_grad.h"
#include <iostream>
#include<Eigen/IterativeLinearSolvers>
#include<Eigen/SparseLU>

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
  
  // Create interpolation matrix
  Eigen::SparseMatrix<double> W, Wx, Wy, Wz, Wx2, Wy2, Wz2, Gx, Gy, Gz;
  // fd_interpolate(nx, ny, nz, h, corner, P, W);

  Eigen::RowVector3d corner1; corner1(0) = corner(0) + h/2;
  fd_interpolate(nx-1, ny, nz, h, corner1, P, Wx);
  Eigen::RowVector3d corner2; corner2(1) = corner(1) + h/2;
  fd_interpolate(nx, ny-1, nz, h, corner2, P, Wy);
  Eigen::RowVector3d corner3; corner3(2) = corner(2) + h/2;
  fd_interpolate(nx, ny, nz-1, h, corner3, P, Wz);

  // Compute total number of grid points
  int ng = nx*ny*nz;

  // Extract total number of given points
  int np = P.rows();

  // Using the triplet insertion method based on Eigen's documentation
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList1, tripletList2, tripletList3;

  // // Extract the x, y, and z parts
  // for (int ii = 0; ii < W.outerSize(); ii++)
  // {
  //     for (Eigen::SparseMatrix<double>::InnerIterator itW(W, ii); itW; ++itW)
  //     {
  //         if (itW.row() < np)
  //           tripletList1.push_back(T(itW.row(), ii, itW.value()));
  //         else if (itW.row() < 2*np)
  //           tripletList2.push_back(T(itW.row()-np, ii, itW.value()));
  //         else
  //           tripletList3.push_back(T(itW.row()-2*np, ii, itW.value()));
  //     }
  // }

  // Extract the x, y, and z parts
  for (int ii = 0; ii < (nx-1)*ny*nz; ii++)
  {
      for (Eigen::SparseMatrix<double>::InnerIterator itW(Wx, ii); itW; ++itW)
      {
          if (itW.row() < np)
            tripletList1.push_back(T(itW.row(), ii, itW.value()));
      }
  }
  for (int ii = 0; ii < nx*(ny-1)*nz; ii++)
  {
      for (Eigen::SparseMatrix<double>::InnerIterator itW(Wy, ii); itW; ++itW)
      {
          if (itW.row() >= np && itW.row() < 2*np)
            tripletList2.push_back(T(itW.row()-np, ii, itW.value()));
      }
  }
  for (int ii = 0; ii < nx*ny*(nz-1); ii++)
  {
      for (Eigen::SparseMatrix<double>::InnerIterator itW(Wz, ii); itW; ++itW)
      {
          if (itW.row() >= 2*np)
            tripletList3.push_back(T(itW.row()-2*np, ii, itW.value()));
      }
  }

  Wx2.resize(np, (nx-1)*ny*nz);
  Wx2.reserve(np*4);
  Wy2.resize(np, nx*(ny-1)*nz);
  Wy2.reserve(np*4);
  Wz2.resize(np, nx*ny*(nz-1));
  Wz2.reserve(np*4);

  Wx2.setFromTriplets(tripletList1.begin(), tripletList1.end());      
  Wy2.setFromTriplets(tripletList2.begin(), tripletList2.end());      
  Wz2.setFromTriplets(tripletList3.begin(), tripletList3.end());      
  
  Wx2.finalize();
  Wy2.finalize();
  Wz2.finalize();

  // Distribute normals on to grid points. First, extract the x, y and z components of the normals
  Eigen::MatrixXd Nx(N.rows(), 1), Ny(N.rows(), 1), Nz(N.rows(), 1);
  for (int ii = 0; ii < N.rows(); ii++)
  {
    Nx(ii, 0) = N(ii, 0);
    Ny(ii, 0) = N(ii, 1);
    Nz(ii, 0) = N(ii, 2);
  }

  // Multiply these vectors with the respective interpolation matrices' transposes
  Eigen::MatrixXd Nvx(ng, 1), Nvy(ng, 1), Nvz(ng, 1);
  Nvx = Wx2.transpose()*Nx;
  Nvy = Wy2.transpose()*Ny;
  Nvz = Wz2.transpose()*Nz;

  // Compute G^TG
  Eigen::SparseMatrix<double> G;
  fd_grad(nx, ny, nz, h, G);

  Eigen::SparseMatrix<double> A = G.transpose()*G;

  // Extract the x, y and z components of the gradient
  // Using the triplet insertion method based on Eigen's documentation
  std::vector<T> tripletList4, tripletList5, tripletList6;

  Gx.resize((nx-1)*ny*nz, nx*ny*nz);
  Gy.resize(nx*(ny-1)*nz, nx*ny*nz);
  Gz.resize(nx*ny*(nz-1), nx*ny*nz);

  // Extract the x, y, and z parts
  for (int ii = 0; ii < G.outerSize(); ii++)
  {
      for (Eigen::SparseMatrix<double>::InnerIterator itG(G, ii); itG; ++itG)
      {
          if (itG.row() < (nx-1)*ny*nz)
            tripletList4.push_back(T(itG.row(), ii, itG.value()));
          else if (itG.row() < (nx-1)*ny*nz + nx*(ny-1)*nz)
            tripletList5.push_back(T(itG.row()-((nx-1)*ny*nz), ii, itG.value()));
          else
            tripletList6.push_back(T(itG.row()-((nx-1)*ny*nz + nx*(ny-1)*nz), ii, itG.value()));
      }
  }

  Gx.setFromTriplets(tripletList4.begin(), tripletList4.end());      
  Gy.setFromTriplets(tripletList5.begin(), tripletList5.end());      
  Gz.setFromTriplets(tripletList6.begin(), tripletList6.end());      
  
  Gx.finalize();
  Gy.finalize();
  Gz.finalize();

  // Now compute the right hand side
  Eigen::VectorXd b = Gx.transpose()*Nvx + Gy.transpose()*Nvy + Gz.transpose()*Nvz;

  // Solve the Ag = b system
  std::cout << "Solving the system..." << std::flush;

  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

  solver.compute(A);
  g = solver.solve(b);

  std::cout << "done." << std::endl;

  // Set iso level - note: this should be done using the formula suggested in the paper; I ran out of time
  double sigma = 0.0;

  // Apply sigma
  for (int ii = 0; ii < nx*ny*nz; ii++)
      g(ii) = g(ii) - sigma;

  // Based on the results, I definitely did not successfully reconstruct the mesh...
  // It looks like I've flipped signs along the way, and / or I think I may have
  // used an incompatible assumption on the ordering of the coordinates (I have
  // assumed they're stored contiguously along x first, then y, then z, but I might have
  // mistakenly been inconsistent). I ran out of time trying to fix this, so I'm
  // submitting it as is (although I definitely want to make it work in the coming days).

  ////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}




