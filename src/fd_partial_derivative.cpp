#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  std::vector< Eigen::Triplet<double> > tripletList;
  std::cout << dir << std::endl;
  int m = D.cols();
  int n = D.rows();
  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      for(int k = 0; k < nz; k++){
        Eigen::RowVector3d coord = Eigen::RowVector3d(i,j,k);
        if(coord(dir) -1 < m && coord(dir) -1 >= 0){
          tripletList.push_back(Eigen::Triplet<double>(i + nx*(j + k * ny), coord(dir)-1, -1));
        }
        if(coord(dir) < m  && coord(dir) >= 0){
          tripletList.push_back(Eigen::Triplet<double>(i + nx*(j + k * ny), coord(dir), 1));
        }
      }
    }
  }
  D.setFromTriplets(tripletList.begin(), tripletList.end());
}
