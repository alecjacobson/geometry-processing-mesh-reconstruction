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
  double px, py, pz, xd, yd, zd, cx, cy, cz;
  int xi, yi, zi;
  for(int i = 0; i < P.rows(); i++) {
    //point coordinates relative to grid-space, in h-units
    px = (P(i,0) - corner(0))/h;
    py = (P(i,1) - corner(1))/h;
    pz = (P(i,2) - corner(2))/h;
    
    //containing voxel's bottom-left-corner coordinates
    cx = std::floor(px);
    cy = std::floor(py);
    cz = std::floor(pz);

    //distances (in h-units)
    xd = px - cx;
    yd = py - cy;
    zd = pz - cz;

    //pushing coeffs for all 8 corners (access other corners by increasing vox-corner index by 1 in each direction)
    triplets.push_back(Eigen::Triplet<double>(i, cx + nx * (cy + ny * cz), (1 - xd) * (1 - yd) * (1 - zd))); //000
    triplets.push_back(Eigen::Triplet<double>(i, (cx + 1) + nx * (cy + ny * cz), xd * (1 - yd) * (1 - zd))); //100
    triplets.push_back(Eigen::Triplet<double>(i, cx + nx * ((cy + 1) + ny * cz), (1 - xd) * yd * (1 - zd))); //010
    triplets.push_back(Eigen::Triplet<double>(i, (cx + 1) + nx * ((cy + 1) + ny * cz), xd * yd * (1 - zd))); //110
    triplets.push_back(Eigen::Triplet<double>(i, cx + nx * (cy + ny * (cz + 1)), (1 - xd) * (1 - yd) * zd)); //001
    triplets.push_back(Eigen::Triplet<double>(i, (cx + 1) + nx * (cy + ny * (cz + 1)), xd * (1 - yd) * zd)); //101
    triplets.push_back(Eigen::Triplet<double>(i, cx + nx * ((cy + 1) + ny * (cz + 1)), (1 - xd) * yd * zd)); //011
    triplets.push_back(Eigen::Triplet<double>(i, (cx + 1) + nx * ((cy + 1) + ny * (cz + 1)), xd * yd * zd)); //111
  }

  W.resize(P.rows(), nx*ny*nz);
  W.setZero();
  W.setFromTriplets(triplets.begin(), triplets.end());

}