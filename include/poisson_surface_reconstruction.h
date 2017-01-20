#ifndef POISSON_SURFACE_RECONSTRUCTION_H
#define POISSON_SURFACE_RECONSTRUCTION_H
#include <Eigen/Core>

// Takes input sample points P and input normals N 
// and gives a watertight mesh using a simplified
// version of [Kazhdan et. al 2006]
//
// Inputs:
//   P  #P by 3 list of input points
//   N  #P by 3 list of input normals associated with each point in P
// Outputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh triangle indinces into V
//
void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F);

#endif
