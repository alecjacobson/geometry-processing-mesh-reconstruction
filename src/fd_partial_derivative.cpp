#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  int nx_c = nx;
  int ny_c = ny;
  int nz_c = nz;
  // Mutable copies

  if (dir == 0) nx_c--;
  else if (dir == 1) ny_c--;
  else nz_c--;

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(nx_c * ny_c * nz_c * 2);

  for (int i = 0; i < nx_c; i++){
    for (int j = 0; j < ny_c; j++){
      for (int k = 0; k < nz_c; k++){
        // In each row
        tripletList.push_back(T(i + j*nx_c + k*nx_c*ny_c, i + j*nx + k*nx*ny, -1));
        
        if (dir == 0){
          tripletList.push_back(T(i + j*nx_c + k*nx_c*ny_c, (i+1) + j*nx + k*nx*ny, 1));
        }
        else if (dir == 1){
          tripletList.push_back(T(i + j*nx_c + k*nx_c*ny_c, i + (j+1)*nx + k*nx*ny, 1));
        }
        else{
          tripletList.push_back(T(i + j*nx_c + k*nx_c*ny_c, i + j*nx + (k+1)*nx*ny, 1));
        }
      }
    }
  }

  D.setFromTriplets(tripletList.begin(), tripletList.end());
}