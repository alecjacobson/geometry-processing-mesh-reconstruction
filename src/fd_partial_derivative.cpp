#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
    //Figure out how dir is used
    
    int dimensions [3];
    
    dimensions[0] = 0;
    dimensions[1] = 0;
    dimensions[2] = 0;
    
    dimensions[dir] = 1;
    
    D.resize((nx-dimensions[0])*(ny-dimensions[1])*(nz-dimensions[2]),nx*ny*nz);
    int pointNo = 0;
    for (int xVal = 0; xVal < nx - dimensions[0] ; xVal ++) {
        for (int yVal = 0; yVal < ny- dimensions[1] ; yVal ++) {
            for (int zVal = 0; zVal < nz - dimensions[2] ; zVal ++) {
                
                //i + nx*(j + k * ny)
                pointNo = xVal + (nx-dimensions[0])*(yVal + zVal*(ny-dimensions[1]));
                
                D.insert(pointNo, xVal + nx*(yVal + zVal*ny)) = - 1.0/h;
                
                
                D.insert(pointNo, (xVal + dimensions[0]) + nx*((yVal + dimensions[1]) + (zVal + dimensions[2])*ny)) = 1.0/h;
                
                
                    
                
                    
        
            }
        }
    
    }
}
