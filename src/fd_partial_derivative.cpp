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
  int x_lim = nx; int y_lim = ny; int z_lim = nz;
  if(dir == 0){
    x_lim -= 1;
  } else if(dir == 1){
    y_lim -= 1;
  } else{
    z_lim -= 1;
  }

  for(int i = 0; i < x_lim; i++){
    for(int j = 0; j < y_lim; j++){
      for(int k = 0; k < z_lim; k++){
        int row_idx = i + x_lim*(j + y_lim*k);
        int col_idx = i + nx*(j + ny*k);
        tripletList.push_back(Eigen::Triplet<double>(row_idx, col_idx, -1));
        int next_col_idx;
        if(dir == 0){
          next_col_idx = (i+1) + nx*(j + ny*k);
        }else if(dir == 1){
          next_col_idx = i + nx*((j+1) + ny*k);
        }else{
          next_col_idx = i + nx*(j + ny*(k+1));
        }
        tripletList.push_back(Eigen::Triplet<double>(row_idx, next_col_idx, 1));
      }
    }
  }
  D.setFromTriplets(tripletList.begin(), tripletList.end());
}
