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
  using Triplet = Eigen::Triplet<double>;

  int Gxs = (nx-1) * (ny  ) * (nz  );
  int Gys = (nx  ) * (ny-1) * (nz  );
  int Gzs = (nx  ) * (ny  ) * (nz-1);

  std::array<int,3> indices = {0,Gxs,Gxs+Gys};
  G.resize(Gxs+Gys+Gzs, nx*ny*nz);

  std::vector<Triplet> trips;
  trips.reserve(2*G.rows());

  for(int i = 0; i < 3; ++i) {
      auto t = fd_partial_derivative_triplets(nx,ny,nz,h,i);
      int ind = indices[i];
      std::transform(t.begin(),t.end(),std::back_inserter(trips), [ind](const Triplet& t) {
              return Triplet{t.row()+ind, t.col(), t.value()};
              });
  }

  G.setFromTriplets(trips.begin(),trips.end());
}
