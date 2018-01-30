#include "fd_partial_derivative.h"
#include <iostream>
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
    //std::cout<<"Here!"<<std::endl;
    if (dir == 0) {
        D.resize((nx-1)*ny*nz, nx*ny*nz);
        for (int i = 0; i < (nx-1); ++i){
            for (int j = 0; j <ny; ++j){
                for (int k = 0; k<nz; ++k){
                    D.insert(i + j * (nx-1) + k * (nx-1)  * ny, i + j * nx + k * nx * ny) = -1;
                    D.insert(i + j * (nx-1) + k * (nx-1)  * ny, i+1 + j * nx + k * nx * ny) = 1;
                }
                    
            }
            
        }
    } else if (dir == 1){
        D.resize(nx*(ny-1)*nz, nx*ny*nz);
        for (int i = 0; i < nx; ++i){
            for (int j = 0; j <(ny -1); ++j){
                for (int k = 0; k<nz; ++k){
                    D.insert(i + j * nx + k * nx  * (ny-1), i + j * nx + k * nx * ny) = -1;
                    D.insert(i + j * nx + k * nx  * (ny-1), i + (j+1) * nx + k * nx * ny) = 1;
                }
                    
            }
            
        }
    } else if (dir == 2){
        D.resize(nx*ny*(nz-1), nx*ny*nz);
        for (int i = 0; i < nx; ++i){
            for (int j = 0; j <ny; ++j){
                for (int k = 0; k<(nz-1); ++k){
                    D.insert(i + j * nx + k * nx  * ny, i + j * nx + k * nx * ny) = -1;
                    D.insert(i + j * nx + k * nx  * ny, i + j * nx + (k+1) * nx * ny) = 1;
                }
                    
            }
            
        }
    }
    //std::cout<<"Done!"<<std::endl;

  ////////////////////////////////////////////////////////////////////////////
}
