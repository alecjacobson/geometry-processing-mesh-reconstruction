#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  int offset[3] = { 0 };
  offset[dir] = 1;
  D.resize((nx-offset[0])*(ny- offset[1])*(nz-offset[2]), nx*ny*nz);
  //From tutorial https://eigen.tuxfamily.org/dox/group__SparseQuickRefPage.html
  //and https://stackoverflow.com/questions/14697375/assigning-a-sparse-matrix-in-eigen
  typedef Eigen::Triplet<double> tuple;
  std::vector<tuple> tupleList;
  tupleList.reserve(2 * D.rows());
  double h_inv = 1 / h;
  for (int i = 0; i < nx-offset[0]; i++) 
    for (int j = 0; j < ny-offset[1]; j++)
      for (int k = 0; k < nz-offset[2]; k++)
      {
        tupleList.push_back(tuple(i + j*(nx- offset[0])+k*(nx- offset[0])*(ny- offset[1]),  i             +  j            *nx +  k            *nx*ny, -h_inv));
        tupleList.push_back(tuple(i + j*(nx- offset[0])+k*(nx- offset[0])*(ny- offset[1]), (i+ offset[0]) + (j+ offset[1])*nx + (k+ offset[2])*nx*ny,  h_inv));
      }
  D.setFromTriplets(tupleList.begin(), tupleList.end());
}
