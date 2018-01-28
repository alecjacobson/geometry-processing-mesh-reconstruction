#include "fd_interpolate.h"

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
  ////////////////////////////////////////////////////////////////////////////
    int n = P.rows();
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(n*8);
    
    for(int i = 0; i < n; i++){
        double coord = P(i,0) - corner(0);
        double mod = coord/h;
        int floorx = int(mod);
        double weightx = 1-(mod - floorx);
        
        coord = P(i,1) - corner(1);
        mod = coord/h;
        int floory = int(mod);
        double weighty = 1-(mod - floory);
        
        coord = P(i,2) - corner(2);
        mod = coord/h;
        int floorz = int(mod);
        double weightz = 1-(mod - floorz);
        
        tripletList.push_back(T(i,floorx + nx*floory + nx*ny*floorz,weightx*weighty*weightz));
        tripletList.push_back(T(i,floorx+1 + nx*floory + nx*ny*floorz,(1-weightx)*weighty*weightz));
        tripletList.push_back(T(i,floorx + nx*(floory+1) + nx*ny*floorz,weightx*(1-weighty)*weightz));
        tripletList.push_back(T(i,floorx+1 + nx*(floory+1) + nx*ny*floorz,(1-weightx)*(1-weighty)*weightz));
        
        tripletList.push_back(T(i,floorx + nx*floory + nx*ny*(floorz+1),weightx*weighty*(1-weightz)));
        tripletList.push_back(T(i,floorx+1 + nx*floory + nx*ny*(floorz+1),(1-weightx)*weighty*(1-weightz)));
        tripletList.push_back(T(i,floorx + nx*(floory+1) + nx*ny*(floorz+1),weightx*(1-weighty)*(1-weightz)));
        tripletList.push_back(T(i,floorx+1 + nx*(floory+1) + nx*ny*(floorz+1),(1-weightx)*(1-weighty)*(1-weightz)));
    }
    W.setFromTriplets(tripletList.begin(), tripletList.end());
}
