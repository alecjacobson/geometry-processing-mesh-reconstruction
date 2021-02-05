#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  double h1=1/h;
  if (dir==0){
    D.resize((nx-1)*ny*nz,nx*ny*nz);
    for (int i=0; i<nx-1; i++){
      //printf("%d ",i);
      for (int j=0; j<ny; j++)
        for (int k=0; k<nz; k++){
          D.insert(i+j*(nx-1)+k*(nx-1)*ny,i+j*nx+k*nx*ny)=-1;
          D.insert(i+j*(nx-1)+k*(nx-1)*ny,i+1+j*nx+k*nx*ny)=1;
        }
    }
  }
  if (dir==1){
    D.resize(nx*(ny-1)*nz,nx*ny*nz);
    for (int i=0; i<nx; i++){
      //printf("%d ",i);
      for (int j=0; j<ny-1; j++)
        for (int k=0; k<nz; k++){
          D.insert(i+j*nx+k*nx*(ny-1),i+j*nx+k*nx*ny)=-1;
          D.insert(i+j*nx+k*nx*(ny-1),i+(j+1)*nx+k*nx*ny)=1;
        }
    }
  }
  if (dir==2){
    D.resize(nx*ny*(nz-1),nx*ny*nz);
    for (int i=0; i<nx; i++){
      //printf("%d ",i);
      for (int j=0; j<ny; j++)
        for (int k=0; k<nz-1; k++){
          D.insert(i+j*nx+k*nx*ny,i+j*nx+k*nx*ny)=-1;
          D.insert(i+j*nx+k*nx*ny,i+j*nx+(k+1)*nx*ny)=1;
        }
    }
  }
}
