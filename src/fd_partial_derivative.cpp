#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  // TODO: refactor...
  if(dir == 0){
    for(int i = 0; i < nx-1; i++){
      for(int j = 0; j < ny; j++){
        for(int k = 0; k < nz; k++){
          for(int col=0; col < nx*ny*nz; col++){
            if(col == i -1){
              D.insert(i + (j*nx) + (k*ny*nx), col) = -1;
            } else if(col == i){
              D.insert(i + (j*nx) + (k*ny*nx), col) = 1;
            }
          }
        }
      }
    }
  }
  if(dir == 1){
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < ny-1; j++){
        for(int k = 0; k < nz; k++){
          for(int col=0; col < nx*ny*nz; col++){
            if(col == j - 1){
              D.insert(i + (j*nx) + (k*ny*nx), col) = -1;
            } else if(col == j){
              D.insert(i + (j*nx) + (k*ny*nx), col) = 1;
            }
          }
        }
      }
    }
  }
  if(dir == 2){
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < ny; j++){
        for(int k = 0; k < nz-1; k++){
          for(int col=0; col < nx*ny*nz; col++){
            if(col == k -1){
              D.insert(i + (j*nx) + (k*ny*nx), col) = -1;
            } else if(col == k){
              D.insert(i + (j*nx) + (k*ny*nx), col) = 1;
            }
          }
        }
      }
    }
  }
}
