#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D) {
    
    // determine the dimensions of the staggered grid. Jump is used to compute the index 
    // of adjacent nodes in the relevant axis based on the derivative direction 
    // (i.e. because of the way the nodes are stored, if we are taking the y-derivative, 
    // we need to jump ahead nx nodes to find the y-neighbour).
    int stx = nx;
    int sty = ny;
    int stz = nz;
    int jump = 1;
    
    // decrement stx, sty, or stz to reflect the number of nodes along the corresponding
    // axis in the staggered grid
    switch (dir) {
        case 0: {
            stx--;
            break;
        }
        
        case 1: {
            sty--;
            jump = nx;
            break;
        }
        
        case 2: {
            stz--;
            jump = nx*ny;
            break;
        }
    }
    
    // resize D so that each row corresponds to a node in the staggered grid and each
    // column corresponds to a node in the regular grid.
    D.resize(stx*sty*stz, nx*ny*nz);
    
    // initialize triplet list used to populate D
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(D.rows()*2);
    
    // loop through each of the nodes in the staggered grid and assign -1 to the regular
    // grid node "left"-neighbour and 1 to the regular grid node "right"-neighbour.
    for (int i = 0; i < stx; i++) {
        for (int j = 0; j < sty; j++) {
            for (int k = 0; k < stz; k++) {
                tripletList.push_back(T(i + stx*(j + sty*k), i + nx*(j + ny*k), -1));
                tripletList.push_back(T(i + stx*(j + sty*k), i + nx*(j + ny*k) + jump, 1)); 
            }
        }
    }
    
    // set D using the triplets list
    D.setFromTriplets(tripletList.begin(), tripletList.end());
}
