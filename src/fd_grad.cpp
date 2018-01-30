#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  std::vector<Eigen::SparseMatrix<double>> D(3);
  for(int i=0; i<D.size(); i++)
    fd_partial_derivative(nx, ny, nz, h, i, D.at(i));

  G.resize(D[0].rows() + D[1].rows() + D[2].rows(), D[0].cols());

  typedef Eigen::Triplet<double> tuple;
  std::vector<tuple> tupleList;
  tupleList.reserve(2*G.rows());

  //Implementation from https://stackoverflow.com/questions/41756428/concatenate-sparse-matrix-eigen
  //for stacking sparse matrix
  int offset[3] = { 0, D[0].rows(), D[0].rows()+D[1].rows()};
  for(int i=0; i<3; i++)
    for (int k = 0; k < D[i].outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(D.at(i), k); it; ++it)
        tupleList.push_back(tuple(offset[i]+it.row(), it.col(), it.value()));

  G.setFromTriplets(tupleList.begin(), tupleList.end());
}

