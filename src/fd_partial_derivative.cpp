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
  int m = 0;
  if (dir == 0) { // nx-1
    int change_x = nx - 1;
    m = change_x*ny*nz;
    D.resize(m, nx*ny*nz);
    for (int i=0; i<change_x;i++){
      for (int j=0; j<ny; j++){
        for (int k=0; k<nz; k++){
          int index = i + nx * (j + ny * k);
          int dindex = i + change_x * (j + ny * k); // change nx
          int nindex = i + 1 + nx * (j + ny * k); // change i
          D.coeffRef(dindex,index) = - 1/h;
          D.coeffRef(dindex,nindex) = 1/h;
        }
      }
    }
  }
  
  else if (dir == 1)
  { // ny - 1
    m = nx*(ny-1)*nz;
    D.resize(m, nx*ny*nz);
    for (int i=0; i<nx;i++){
      for (int j=0; j<ny-1; j++){
        for (int k=0; k<nz; k++){
          int index = i + nx * (j + ny * k);
          int dindex = i + nx * (j + (ny-1) * k); // change ny
          int nindex = i + nx * (j + 1 + ny * k); // change j
          D.coeffRef(dindex,index) = - 1/h;
          D.coeffRef(dindex,nindex) = 1/h;
        }
      }
    }
  }
  
  else { // dir == 2 nz-1
    m = nx*ny*(nz-1);
    D.resize(m, nx*ny*nz);
    for (int i=0; i<nx;i++){
      for (int j=0; j<ny; j++){
        for (int k=0; k<nz-1; k++){
          int index = i + nx * (j + ny * k);
          int dindex = i + nx * (j + ny * k); // change nz
          int nindex = i + nx * (j + ny * (k+1)); // change z
          D.coeffRef(dindex,index) = - 1/h;
          D.coeffRef(dindex,nindex) = 1/h;
        }
      }
    }
  }

 }
