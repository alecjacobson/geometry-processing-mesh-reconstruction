#include "fd_interpolate.h"
#include <iostream>

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
    //NEED TO FIGURE OUT INDEXING
  ////////////////////////////////////////////////////////////////////////////
    
    
    int dims[3];
    dims[0] = nx;
    dims[1] = ny;
    dims[2] = nz;
    W.resize(P.rows(), nx*ny*nz);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(8*P.rows());
    
    
    double P_off[3], P_mod[3];
    for (int pointNo = 0; pointNo < P.rows(); pointNo++) {
        for (int P_ind = 0; P_ind < 3; P_ind ++) {
            
            P_mod[P_ind] = floor((double) (P(pointNo,P_ind) - corner(P_ind)) / h);
            P_off[P_ind] = fmod(P(pointNo,P_ind) - corner(P_ind), h);
            
        }
        
        //Do the indexing
        int curIndices[3], ind_val, tempVal;
        
        for (int curInd = 0; curInd < 8; curInd++) {
            double curWeight = 1;
            ind_val = curInd;
            for (int expVal = 2; expVal > -1; expVal --) {
                tempVal  = floor((double) ind_val / pow(2, expVal));
                curIndices[expVal] = P_mod[expVal] + tempVal ;
                
                if (tempVal == 1) {
                    curWeight = curWeight * (1 - P_off[expVal]);
                    
                } else {
                    curWeight = curWeight * P_off[expVal];
                
                
                }
                ind_val = ind_val - tempVal * pow(2,2);
                
                 //std::cout << "curInd: " << curInd << " Offset" << P_off[expVal] << " " << curWeight << "\n";
                
            }
            //i + nx*(j + k * ny)
           
            
            tripletList.push_back(T(pointNo, curIndices[0] + nx*(curIndices[1] + curIndices[2]*ny),curWeight));
            
        }
    
        
    }

    W.setFromTriplets(tripletList.begin(), tripletList.end());
    
}

