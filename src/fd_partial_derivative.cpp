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
  ////////////////////////////////////////////////////////////////////////////
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    
    int m = nx*ny*nz;

    for(int l = 0; l < m+1; l++){
        int loc1 = l;
        int loc2 = l+1;
        int limit = (nx-1)*ny*nz;
        if (dir == 1) {
            loc2 = nx+(l);
            limit = nx*(ny-1)*nz;
        } else if (dir == 2){
            loc2 = nx*ny+(l);
            limit = nx*ny*(nz-1);
        }
        if(loc1<limit)
            tripletList.push_back(T(loc1,l,1));
        if(loc2<limit)
            tripletList.push_back(T(loc2,l,-1));
    }
    D.setFromTriplets(tripletList.begin(), tripletList.end());
}
