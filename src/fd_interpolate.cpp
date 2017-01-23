#include "fd_interpolate.h"
#include <math.h>
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
  // set size of W
  W.resize(P.rows(), nx*ny*nz);
  W.reserve(P.rows() * 8);
  
  double xd, yd, zd, ci, cj, ck;
  int a, b, c, ind;
  
  for(int p = 0; p < P.rows(); p++) {
    // steps to c000
    a = (P(p,0) - corner(0,0)) / h;
    b = (P(p,1) - corner(0,1)) / h;
    c = (P(p,2) - corner(0,2)) / h;
    
    // find c000 value
    ci = a * h + corner(0,0);
    cj = b * h + corner(0,1);
    ck = c * h + corner(0,2);
    
    // find xd, yd and zd;
    xd = (P(p,0) - ci) / h;
    yd = (P(p,1) - cj) / h;
    zd = (P(p,2) - ck) / h;
    
    // c000 weight
    ind = a + nx*(b + c * ny);
    W.insert(p, ind) = (1 - xd)*(1 - yd)*(1 - zd);
    
    // c100 weight
    ind = (a + 1) + nx*(b + c * ny);
    W.insert(p, ind) = xd*(1 - yd)*(1 - zd);
    
    // c010 weight
    ind = a + nx*( (b + 1) + c * ny);
    W.insert(p, ind) = (1-xd)*yd*(1-zd);
    
    // c110 weight
    ind = (a + 1) + nx*( (b + 1) + c * ny);
    W.insert(p, ind) = xd * yd * (1-zd);
    
    // c001 weight
    ind = a + nx*(b + (c + 1) * ny);
    W.insert(p, ind) = (1-xd)*(1-yd)*zd;
    
    // c101 weight
    ind = (a + 1) + nx*(b + (c + 1) * ny);
    W.insert(p, ind) = xd * (1-yd) * zd;
    
    // c011 weight
    ind = a + nx*( (b + 1) + (c + 1) * ny);
    W.insert(p, ind) = (1 - xd) * yd * zd;
    
    // c111 weight
    ind = (a + 1) + nx*((b + 1) + (c + 1) * ny);
    W.insert(p, ind) = xd * yd * zd;
  }
}
