#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <iostream>

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
    
    //Use Reserve
    G.resize((nx-1)*ny*nz+ nx*(ny-1)*nz+ nx*ny*(nz-1),nx*ny*nz);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(2*((nx-1)*ny*nz+ nx*(ny-1)*nz+ nx*ny*(nz-1)));
    int dims[3];
    dims[0] = nx;
    dims[1] = ny;
    dims[2] = nz;
    
    int counter = 0, curRow, curCol, curVal;
    for (int dir = 0; dir < 3; dir ++) {
        Eigen::SparseMatrix<double> tempMat;
        std::cout << "started Grad " << dir << "\n";
        fd_partial_derivative(nx,ny,nz,h,dir, tempMat);
        std::cout << "Grad " << dir << "\n";
        //Temporary solution cite eigen page
        for (int k=0; k<tempMat.outerSize(); ++k){
            for (Eigen::SparseMatrix<double>::InnerIterator it(tempMat,k); it; ++it)
            {
                curVal = it.value();
                curRow = it.row();   // row index
                curCol = it.col();   // col index (here it is equal to k)
                
                tripletList.push_back(T(curRow + counter,curCol,curVal));
                
                //std::cout << curRow + counter << "\n";
            }
            //std::cout << k << "\n";
    
        }
        
        //incrementing count
        dims[dir] = dims[dir] - 1;
        counter = counter + dims[0]*dims[1]*dims[2];
        dims[dir] = dims[dir] + 1;
        
        
    }
    G.setFromTriplets(tripletList.begin(), tripletList.end());
    //std::cout << (nx-1)*ny*nz+ nx*(ny-1)*nz+ nx*ny*(nz-1) << "\n";
    //std::cout << G.nonZeros() << "\n";
}
