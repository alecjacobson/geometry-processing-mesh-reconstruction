#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include "fd_interpolate.h"
#include "fd_grad.h"
#include "fd_partial_derivative.h"
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
  Eigen::VectorXd g = Eigen::VectorXd::Ones(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
    
    
    const int vDim =(nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1);
    Eigen::MatrixXd smallV;
    smallV.resize(vDim,1);
    Eigen::SparseMatrix<double> G;
        
    double sigma = 0;
    int dims[3];
    dims[0] = nx;
    dims[1] = ny;
    dims[2] = nz;
    int counter = 0;
    
    //Compute v_x, v_y, v_z
    for (int dir = 0; dir < 3; dir ++) {
        Eigen::SparseMatrix<double> tempMat;
        Eigen::MatrixXd tempV;
        
        //Shifts the box for the appropriate updated dimensions
        dims[dir] = dims[dir] - 1;
        corner(dir) = corner(dir) + h;
        fd_interpolate(dims[0], dims[1],dims[2],h,corner,P, tempMat);
        corner(dir) = corner(dir) - h;

        tempV = tempMat.transpose()* N.col(dir);
        

        
        for (int xVal = 0; xVal < dims[0] ; xVal ++) {
            for (int yVal = 0; yVal < dims[1]; yVal ++) {
                for (int zVal = 0; zVal <  dims[2] ; zVal ++) {
                    
                    int pointNo = xVal + (dims[0])*(yVal + zVal*(dims[1]));
                    
                    //Computing the values for v
                    
                    
                    
                    smallV(counter + pointNo, 0) = tempV(pointNo,0);
                    
                }
            }
        }
        
        
        counter = counter + dims[0]*dims[1]*dims[2];
        dims[dir] = dims[dir] + 1;
        
        
    }
   
    //Compute the matrix of concatenated partial derivatives
    fd_grad(nx,ny,nz,h,G);
    
    //Use Conjugate Gradient solver
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    
    Eigen::SparseMatrix<double> tempG = G.transpose()*G;
    
    
    solver.compute(tempG);
    
    //Solves GtGg = G^tv for g;
    g = solver.solve(G.transpose()*smallV);

    
    Eigen::SparseMatrix<double> tempMat;
    fd_interpolate(dims[0], dims[1],dims[2],h,corner,P, tempMat);
    
    Eigen::MatrixXd newG = tempMat * g;
    for (int i = 0; i < P.rows(); i ++){
    sigma +=  newG(i,0);
    }
    
    sigma = sigma / (double) P.rows();
    g.array() -= sigma;
    
  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
