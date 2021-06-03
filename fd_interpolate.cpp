#include "fd_interpolate.h"

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W) {

    // resize W so that each row corresponds to a point in P and each column corresponds
    // to a node in the grid. Each row of W will have 8 non-zero elements (the weights
    // for each of the 8 nodes comprising the voxel containing the given point).
    W.resize(P.rows(), nx*ny*nz);
    
    // initialize a triplet list used to populate W
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(W.rows()*8);
    
    // compute the weights for each point in the point cloud
    for (int p = 0; p < P.rows(); p++) {
        
        // convert the cartesian coordinates of the point to coordinates relative to
        // the grid front-bottom-left corner in units of h
        double x = (P(p,0) - corner(0))/h;
        double y = (P(p,1) - corner(1))/h;
        double z = (P(p,2) - corner(2))/h;
        
        // h is the distance between consecutive grid nodes along a given
        // axis, therefore we can use the floor operator to determine the "axis index"
        // of the front-bottom-left node of the voxel in which the point is situated.
        int vx = floor(x);
        int vy = floor(y);
        int vz = floor(z);
        
        // compute the remainder (i.e. the distance in each axis from the point to the 
        // front-bottom-left node).
        double rx = x - vx;
        double ry = y - vy;
        double rz = z - vz;
        
        // compute the weights for the nodes comprising the relevant voxel and add them
        // to the triplet list
        tripletList.push_back(T(p, vx + nx*vy + nx*ny*vz, (1-rx)*(1-ry)*(1-rz)));
        tripletList.push_back(T(p, vx + nx*vy + nx*ny*(vz + 1), (1-rx)*(1-ry)*(rz)));
        tripletList.push_back(T(p, vx + nx*(vy + 1) + nx*ny*vz, (1-rx)*(ry)*(1-rz)));
        tripletList.push_back(T(p, vx + nx*(vy + 1) + nx*ny*(vz + 1), (1-rx)*(ry)*(rz)));
        tripletList.push_back(T(p, (vx + 1) + nx*vy + nx*ny*vz, (rx)*(1-ry)*(1-rz)));
        tripletList.push_back(T(p, (vx + 1) + nx*vy + nx*ny*(vz + 1), (rx)*(1-ry)*(rz)));
        tripletList.push_back(T(p, (vx + 1) + nx*(vy + 1) + nx*ny*vz, (rx)*(ry)*(1-rz)));
        tripletList.push_back(T(p, (vx + 1) + nx*(vy + 1) + nx*ny*(vz + 1), rx*ry*rz));
                                        
    }
    
    // populate W using the triplets list
    W.setFromTriplets(tripletList.begin(), tripletList.end());
}
