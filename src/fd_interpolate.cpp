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
  std::vector<Eigen::Triplet<double>> triplets; 
  double px, py, pz, xd, yd, zd, cx, cy, cz; // split everything up to think it through properly...
  int xi, yi, zi;
  for(int i = 0; i < P.rows(); i++) {
    //current  position relative to grid corner
    px = P(i,0) - corner(0);
    py = P(i,1) - corner(1);
    pz = P(i,2) - corner(2);

    //indices in grid-space
    xi = std::floor(px/h);
    yi = std::floor(py/h);
    zi = std::floor(pz/h);
    
    //position of voxel-corner at that index
    cx = xi*h;
    cy = yi*h;
    cz = zi*h;

    //distances from point to voxel-corner, to calculate weights
    xd = px - cx;
    yd = py - cy;
    zd = pz - cz;

    //pushing coeffs for all 8 corners (access other corners by increasing vox-corner index by 1 in each direction)
    triplets.push_back(Eigen::Triplet<double>(i, xi + nx * (yi + ny * zi), (1 - xd) * (1 - yd) * (1 - zd))); //000
    triplets.push_back(Eigen::Triplet<double>(i, (xi + 1) + nx * (yi + ny * zi), xd * (1 - yd) * (1 - zd))); //100
    triplets.push_back(Eigen::Triplet<double>(i, xi + nx * ((yi + 1) + ny * zi), (1 - xd) * yd * (1 - zd))); //010
    triplets.push_back(Eigen::Triplet<double>(i, (xi + 1) + nx * ((yi + 1) + ny * zi), xd * yd * (1 - zd))); //110
    triplets.push_back(Eigen::Triplet<double>(i, xi + nx * (yi + ny * (zi + 1)), (1 - xd) * (1 - yd) * zd)); //001
    triplets.push_back(Eigen::Triplet<double>(i, (xi + 1) + nx * (yi + ny * (zi + 1)), xd * (1 - yd) * zd)); //101
    triplets.push_back(Eigen::Triplet<double>(i, xi + nx * ((yi + 1) + ny * (zi + 1)), (1 - xd) * yd * zd)); //011
    triplets.push_back(Eigen::Triplet<double>(i, (xi + 1) + nx * ((yi + 1) + ny * (zi + 1)), xd * yd * zd)); //111
  }

  W.resize(P.rows(), nx*ny*nz);
  W.setZero();
  W.setFromTriplets(triplets.begin(), triplets.end());

}