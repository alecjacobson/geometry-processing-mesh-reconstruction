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

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(8*P.rows());

  // loop for each point
  for (int l = 0; l < P.rows(); ++l) {
    const double& px = P(l,0);
    const double& py = P(l,1);
    const double& pz = P(l,2);  

    int i0 = (px - corner(0))/h;
    int j0 = (py - corner(1))/h;
    int k0 = (pz - corner(2))/h;
    int i1 = i0 + 1;
    int j1 = j0 + 1;
    int k1 = k0 + 1;

    double xd = (px - i0*h - corner(0))/h;
    double yd = (py - j0*h - corner(1))/h;
    double zd = (pz - k0*h - corner(2))/h;

    double val = (1.0 - xd)*(1.0 - yd)*(1.0 - zd);
    tripletList.push_back(T(l, i0 + j0*nx + k0*ny*nx, val));

    val = xd*(1.0 - yd)*(1.0 - zd);
    tripletList.push_back(T(l, i1 + j0*nx + k0*ny*nx, val));

    val = (1.0 - xd)*yd*(1.0-zd);
    tripletList.push_back(T(l, i0 + j1*nx + k0*ny*nx, val));

    val = xd*yd*(1.0 - zd);
    tripletList.push_back(T(l, i1 + j1*nx + k0*ny*nx, val));

    val = (1.0 - xd)*(1.0 - yd)*zd;
    tripletList.push_back(T(l, i0 + j0*nx + k1*ny*nx, val));

    val = xd*(1.0 - yd)*zd;
    tripletList.push_back(T(l, i1 + j0*nx + k1*ny*nx, val));

    val = (1.0 - xd)*yd*zd;
    tripletList.push_back(T(l, i0 + j1*nx + k1*ny*nx, val));

    val = xd*yd*zd;
    tripletList.push_back(T(l, i1 + j1*nx + k1*ny*nx, val));

  } // end loop l

  W = Eigen::SparseMatrix<double>(P.rows(), nx*ny*nz);
  W.setFromTriplets(tripletList.begin(), tripletList.end());
}
