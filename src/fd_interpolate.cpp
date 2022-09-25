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
    typedef Eigen::Triplet<double> Tr;
    std::vector<Tr> trips;
    for (int i = 0; i < P.rows(); i++) {
        Eigen::RowVector3d p = (P.row(i) - corner) / h;
        Eigen::RowVector3d fl;
        fl << std::floor(p[0]), std::floor(p[1]), std::floor(p[2]);
        Eigen::RowVector3d d = p - fl;
        trips.push_back(Tr(i, fl[0] + fl[1]*nx + fl[2]*nx*ny, (1-d[0])*(1-d[1])*(1-d[2])));
        trips.push_back(Tr(i, (fl[0]+1) + fl[1]*nx + fl[2]*nx*ny, d[0]*(1-d[1])*(1-d[2])));
        trips.push_back(Tr(i, fl[0] + (fl[1]+1)*nx + fl[2]*nx*ny, (1-d[0])*d[1]*(1-d[2])));
        trips.push_back(Tr(i, fl[0] + fl[1]*nx + (fl[2]+1)*nx*ny, (1-d[0])*(1-d[1])*d[2]));
        trips.push_back(Tr(i, (fl[0]+1) + (fl[1]+1)*nx + fl[2]*nx*ny, d[0]*d[1]*(1-d[2])));
        trips.push_back(Tr(i, (fl[0]+1) + fl[1]*nx + (fl[2]+1)*nx*ny, d[0]*(1-d[1])*d[2]));
        trips.push_back(Tr(i, fl[0] + (fl[1]+1)*nx + (fl[2]+1)*nx*ny, (1-d[0])*d[1]*d[2]));
        trips.push_back(Tr(i, (fl[0]+1) + (fl[1]+1)*nx + (fl[2]+1)*nx*ny, d[0]*d[1]*d[2]));
    }
    W.setFromTriplets(trips.begin(), trips.end());
}
