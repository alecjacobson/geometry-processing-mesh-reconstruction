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

    W.resize(P.rows(), nx * ny * nz);

    std::vector<Eigen::Triplet<double>> buffer;

    for(int i = 0; i < P.rows(); i++){
        int x0 = (P(i, 0) - corner(0)) / h;
        int x1 = x0 + 1;
        double dx = (P(i, 0) - (x0 * h + corner(0))) / h;

        int y0 = (P(i, 1) - corner(1)) / h;
        int y1 = y0 + 1;
        double dy = (P(i, 1) - (y0 * h + corner(1))) / h;

        int z0 = (P(i, 2) - corner(2)) / h;
        int z1 = z0 + 1;
        double dz = (P(i, 2) - (z0 * h + corner(2))) / h;

        buffer.push_back(Eigen::Triplet<double>(i, x0 + y0*nx + z0*ny*nx,(1.0 - dx)*(1.0 - dy)*(1.0 - dz)));
        buffer.push_back(Eigen::Triplet<double>(i, x1 + y0*nx + z0*ny*nx,(      dx)*(1.0 - dy)*(1.0 - dz)));
        buffer.push_back(Eigen::Triplet<double>(i, x0 + y1*nx + z0*ny*nx,(1.0 - dx)*(      dy)*(1.0 - dz)));
        buffer.push_back(Eigen::Triplet<double>(i, x1 + y1*nx + z0*ny*nx,(      dx)*(      dy)*(1.0 - dz)));
        buffer.push_back(Eigen::Triplet<double>(i, x0 + y0*nx + z1*ny*nx,(1.0 - dx)*(1.0 - dy)*(      dz)));
        buffer.push_back(Eigen::Triplet<double>(i, x1 + y0*nx + z1*ny*nx,(      dx)*(1.0 - dy)*(      dz)));
        buffer.push_back(Eigen::Triplet<double>(i, x0 + y1*nx + z1*ny*nx,(1.0 - dx)*(      dy)*(      dz)));
        buffer.push_back(Eigen::Triplet<double>(i, x1 + y1*nx + z1*ny*nx,(      dx)*(      dy)*(      dz)));
    }

    W.setFromTriplets(buffer.begin(), buffer.end());
}
