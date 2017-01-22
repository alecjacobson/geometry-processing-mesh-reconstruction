#include "fd_interpolate.h"
#include <vector>

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

    using Triplet = Eigen::Triplet<double>;

    auto indexer = [&](int i, int j, int k) -> int {
        return i + nx*(j + k * ny);
    };
    auto indexer_v = [&](const Eigen::RowVectorXi& ijk, int i=0, int j=0, int k=0) -> int {
        return indexer(ijk.x()+i,ijk.y()+j,ijk.z()+k);
    };


    std::vector<Triplet> trips;
    trips.reserve(8*P.cols());


    for(int r = 0; r < P.rows(); ++r) {
        auto pr = P.row(r);

        //material space coordinates
        Eigen::RowVector3d mc = (pr - corner) / h;

        //index values
        //TODO: figure out why ArrayBase::floor doesn't exist?
        Eigen::RowVector3i iv = mc.unaryExpr([](double d){return floor(d);}).cast<int>();
        iv(0) = std::min(nx-1,np.max(0,iv(0)));
        iv(1) = std::min(ny-1,np.max(0,iv(1)));
        iv(2) = std::min(nz-1,np.max(0,iv(2)));

        //Interpolation alphas
        Eigen::Matrix<double,2,3> a;
        a.row(1) = mc - iv.cast<double>();
        a.row(0) = 1 - a.row(1).array();




        for(int i = 0; i < 2; ++i) {
            for(int j = 0; j < 2; ++j) {
                for(int k = 0; k < 2; ++k) {

                    int ind = indexer_v(iv.cast<int>(),i,j,k);
                    trips.emplace_back(
                            Triplet{r,ind,a(i,0)*a(j,1)*a(k,2)
                            });
                }
            }
        }


    }
    W.resize(P.rows(),nx*ny*nz);

    W.setFromTriplets(trips.begin(),trips.end());


}
