#include "fd_interpolate.h"
#include <iostream>

using namespace std;

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
  int numPt = P.rows();
  std::vector<Eigen::Triplet<double> > triplets(8*numPt);

  Eigen::RowVector3d cornerToPt;
  double voxelVec_x, voxelVec_y, voxelVec_z;
  int xIdx, yIdx, zIdx;
  int k = 0;
  for (int ii = 0; ii < numPt; ii++){
    cornerToPt = P.row(ii) - corner;
    // cout << cornerToPt << endl;
     
    // W index
    xIdx = int(floor(cornerToPt(0) / h));
    yIdx = int(floor(cornerToPt(1) / h));
    zIdx = int(floor(cornerToPt(2) / h));
    // cout << xIdx << yIdx << zIdx << endl;

    // normalize vec
    voxelVec_x = fmod(cornerToPt(0),h) / h;
    voxelVec_y = fmod(cornerToPt(1),h) / h;
    voxelVec_z = fmod(cornerToPt(2),h) / h;
    // cout << voxelVec_x << voxelVec_y << voxelVec_z << endl;

    // W weight
    // W.coeffRef(ii, xIdx + yIdx*nx + zIdx*ny*nx) = (1-voxelVec_x) * (1-voxelVec_y) * (1-voxelVec_z); // W0
    // W.coeffRef(ii, (xIdx+1) + yIdx*nx + zIdx*ny*nx) = voxelVec_x * (1-voxelVec_y) * (1-voxelVec_z); // W1
    // W.coeffRef(ii, xIdx + (yIdx+1)*nx + zIdx*ny*nx) = (1-voxelVec_x) * voxelVec_y * (1-voxelVec_z); // W2
    // W.coeffRef(ii, (xIdx+1) + (yIdx+1)*nx + zIdx*ny*nx) = voxelVec_x * voxelVec_y * (1-voxelVec_z); // W3
    // W.coeffRef(ii, xIdx + yIdx*nx + (zIdx+1)*ny*nx) = (1-voxelVec_x) * (1-voxelVec_y) * voxelVec_z; // W4
    // W.coeffRef(ii, (xIdx+1) + yIdx*nx + (zIdx+1)*ny*nx) = voxelVec_x * (1-voxelVec_y) * voxelVec_z; // W5
    // W.coeffRef(ii, xIdx + (yIdx+1)*nx + (zIdx+1)*ny*nx) = (1-voxelVec_x) * voxelVec_y * voxelVec_z; // W6
    // W.coeffRef(ii, (xIdx+1) + (yIdx+1)*nx + (zIdx+1)*ny*nx) = voxelVec_x * voxelVec_y * voxelVec_z; // W7

    // construct triplet
    triplets[k++] = Eigen::Triplet<double>(ii, xIdx + yIdx*nx + zIdx*ny*nx, (1-voxelVec_x) * (1-voxelVec_y) * (1-voxelVec_z));
    triplets[k++] = Eigen::Triplet<double>(ii, (xIdx+1) + yIdx*nx + zIdx*ny*nx, voxelVec_x * (1-voxelVec_y) * (1-voxelVec_z));
    triplets[k++] = Eigen::Triplet<double>(ii, xIdx + (yIdx+1)*nx + zIdx*ny*nx, (1-voxelVec_x) * voxelVec_y * (1-voxelVec_z));
    triplets[k++] = Eigen::Triplet<double>(ii, (xIdx+1) + (yIdx+1)*nx + zIdx*ny*nx, voxelVec_x * voxelVec_y * (1-voxelVec_z));
    triplets[k++] = Eigen::Triplet<double>(ii, xIdx + yIdx*nx + (zIdx+1)*ny*nx, (1-voxelVec_x) * (1-voxelVec_y) * voxelVec_z);
    triplets[k++] = Eigen::Triplet<double>(ii, (xIdx+1) + yIdx*nx + (zIdx+1)*ny*nx, voxelVec_x * (1-voxelVec_y) * voxelVec_z);
    triplets[k++] = Eigen::Triplet<double>(ii, xIdx + (yIdx+1)*nx + (zIdx+1)*ny*nx, (1-voxelVec_x) * voxelVec_y * voxelVec_z);
    triplets[k++] = Eigen::Triplet<double>(ii, (xIdx+1) + (yIdx+1)*nx + (zIdx+1)*ny*nx, voxelVec_x * voxelVec_y * voxelVec_z);
  }
  W.resize(P.rows(),nx*ny*nz);
  W.setFromTriplets(triplets.begin(), triplets.end());
}
