#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
    Eigen::SparseMatrix<double> tempMat;
    
    int dims[3];
    dims[0] = nx;
    dims[1] = ny;
    dims[2] = nz;
    
    int counter = 0, curRow, curCol, curVal;
    for (int dir = 0; dir < 3; dir ++) {
        fd_partial_derivative(nx,ny,nz,h,dir, tempMat);
        
        //Temporary solution cite eigen page
        for (int k=0; k<tempMat.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(tempMat,k); it; ++it)
            {
                curVal = it.value();
                curRow = it.row();   // row index
                curCol = it.col();   // col index (here it is equal to k)
                
                
                G.insert(curRow + counter, curCol) = curVal;
            }
    
        }
        //incrementing count
        dims[dir] = dims[dir] - 1;
        counter = counter + dims[0]*dims[1]*dims[2];
        dims[dir] = dims[dir] + 1;
        
        
    }
}
