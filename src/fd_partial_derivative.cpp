#include "fd_partial_derivative.h"

int calIdx(
  const int xIdx, const int yIdx, const int zIdx,
  const int nx, const int ny, const int nz){
  return xIdx + nx*yIdx + zIdx*nx*ny;
}

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
  double colnx = nx;
  double colny = ny;
  double colnz = nz;
  if (dir == 0)
    colnx -= 1;
  else if (dir == 1)
    colny -= 1;
  else
    colnz -= 1;

  D.resize(colnx*colny*colnz, nx*ny*nz);
  int colIdx, rowIdx, nextIdx;
  for(int ii = 0; ii < colnx; ii++){
    for(int jj = 0; jj < colny; jj++){
      for(int kk = 0; kk < colnz; kk++){
        rowIdx = calIdx(ii, jj, kk, colnx, colny, colnz);
        colIdx = calIdx(ii, jj, kk, nx, ny, nz);
        if (dir == 0){
          D.coeffRef(rowIdx, colIdx) = -1/h;
          nextIdx = calIdx(ii+1, jj, kk, nx, ny, nz);
          D.coeffRef(rowIdx, nextIdx) = 1/h;
        }
        else if (dir == 1){
          D.coeffRef(rowIdx, colIdx) = -1/h;
          nextIdx = calIdx(ii, jj+1, kk, nx, ny, nz);
          D.coeffRef(rowIdx, nextIdx) = 1/h;
        }
        else{
          D.coeffRef(rowIdx, colIdx) = -1/h;
          nextIdx = calIdx(ii, jj, kk+1, nx, ny, nz);
          D.coeffRef(rowIdx, nextIdx) = 1/h;
        }
      }
    }
  }
}
