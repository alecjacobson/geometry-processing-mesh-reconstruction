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
    smallV = Eigen::MatrixXd::Ones(vDim,1);
    Eigen::SparseMatrix<double> G;
    Eigen::MatrixXd OnesMat;
    OnesMat = Eigen::MatrixXd::Ones(1,P.rows());
    
    double sigma = 0;
    int dims[3];
    dims[0] = nx;
    dims[1] = ny;
    dims[2] = nz;
    int counter = 0;
    
    for (int dir = 0; dir < 3; dir ++) {
        Eigen::SparseMatrix<double> tempMat;
        Eigen::MatrixXd tempV;
        dims[dir] = dims[dir] - 1;
        
        fd_interpolate(dims[0], dims[1],dims[2],h,corner,P, tempMat);
        
        tempV.resize(dims[0]*dims[1]*dims[2],1);
        tempV = tempMat.transpose()* N.col(dir);
        
        /*for (int i = 0; i < N.rows(); i ++){
            std::cout << N(i,dir) << "\n";
        }
        
        
        for (int k=0; k<tempMat.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(tempMat,k); it; ++it)
            {
                double curVal = it.value();
                std::cout << curVal << "\n";
            }
            //std::cout << k << "\n";
            
        }*/
        
        
        for (int i = 0; i < dims[0]*dims[1]*dims[2]; i ++){
            
            smallV(counter + i, 0) = tempV(i,0);
        
        }
        
        counter = counter + dims[0]*dims[1]*dims[2];
        dims[dir] = dims[dir] + 1;
        
        
    }
   
    
    fd_grad(nx,ny,nz,h,G);
    
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    
    Eigen::SparseMatrix<double> tempG = G.transpose()*G;
    
    //Check tempG values
    double curVal,curRow, curCol;
    for (int k=0; k<smallV.rows();k++) {
        std::cout << smallV(k,0) << "\n";
    
    }
    //Delete above after
    
    solver.compute(tempG);
    std::cout << tempG.nonZeros() << "hi\n";
    std::cout << nx << " "<< ny << " "<< nz << " "<< "\n";
    g = solver.solve(G.transpose()*smallV);
    std::cout << "hi2\n";
    
    Eigen::SparseMatrix<double> tempMat;
    fd_interpolate(dims[0], dims[1],dims[2],h,corner,P, tempMat);
    
    Eigen::MatrixXd newG = tempMat * g;
    for (int i = 0; i < P.rows(); i ++){
    sigma +=  newG(i,0);
    }
    
    sigma = sigma / (double) P.rows();
    g.array() -= sigma;
    
    /*for (int i = 0; i < g.rows() ;i ++){
        std::cout << g(i) << "\n";
    }*/
  //Compute W Matrices here.
    
    
  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
