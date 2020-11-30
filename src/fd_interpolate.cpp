#include "fd_interpolate.h"
#include <math.h>
#include <vector>

#define IND(i,j,k) i + nx*(j + (k) * ny)

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  W.resize(P.rows(), nx*ny*nz);
  
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  
  for (int i = 0; i < W.rows(); i++) {
  
    double px = P(i,0);
    double py = P(i,1);
    double pz = P(i,2);
    
    // obtain corresponding (non-integer) grid coordinates
    double x = (px - corner[0]) / h;
    double y = (py - corner[1]) / h;
    double z = (pz - corner[2]) / h;
    
    // obtain the coordinates involved in points of surrounding cube
    int x_l = floor(x);
    int x_r = ceil(x);
    int y_l = floor(y);
    int y_r = ceil(y);
    int z_l = floor(z);
    int z_r = ceil(z);
    
    double x_d = x - x_l;
    double y_d = y - y_l;
    double z_d = z - z_l;
    
    // set the weights for each of the eight surrounding points
    tripletList.push_back(T(i, IND(x_l,y_l,z_l), (1 - x_d) * (1 - y_d) * (1 - z_d)));
    tripletList.push_back(T(i, IND(x_l,y_l,z_r), (1 - x_d) * (1 - y_d) * z_d));
    tripletList.push_back(T(i, IND(x_l,y_r,z_l), (1 - x_d) * y_d * (1 - z_d)));
    tripletList.push_back(T(i, IND(x_l,y_r,z_r), (1 - x_d) * y_d * z_d));
    tripletList.push_back(T(i, IND(x_r,y_l,z_l), x_d * (1 - y_d) * (1 - z_d)));
    tripletList.push_back(T(i, IND(x_r,y_l,z_r),  x_d * (1 - y_d) * z_d));
    tripletList.push_back(T(i, IND(x_r,y_r,z_l), x_d * y_d * (1 - z_d)));
    tripletList.push_back(T(i, IND(x_r,y_r,z_r), x_d * y_d * z_d));
  }
  
  W.setFromTriplets(tripletList.begin(), tripletList.end());
}
