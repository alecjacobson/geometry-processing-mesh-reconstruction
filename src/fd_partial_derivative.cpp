#include "fd_partial_derivative.h"

class Index
{
  public:
    int nx, ny, nz;
    Index(int x, int y, int z){
        nx = x; ny = y; nz = z;
    };
    int getIndex(int i, int j, int k){
      return i + nx*j + nx*ny*k;
    };
};

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

  Index indexer(nx_c, ny_c, nz_c);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(nx_c * ny_c * nz_c * 2);

  for (int i = 0; i < nx_c; i++){
    for (int j = 0; j < ny_c; j++){
      for (int k = 0; k < nz_c; k++){
        // In each row
        tripletList.push_back(T(indexer.getIndex(i, j, k), indexer.getIndex(i, j, k), -1./h));
        // tripletList.push_back(T(i + j*nx + k*nx*ny, i + j*nx + k*nx*ny, -1./h));
        
        if (dir == 0){
          tripletList.push_back(T(indexer.getIndex(i, j, k), indexer.getIndex(i+1, j, k), 1./h));
          // tripletList.push_back(T(i + j*nx + k*nx*ny, (i+1) + j+*nx + k*nx*ny, 1./h));
        }
        else if (dir == 1){
          tripletList.push_back(T(indexer.getIndex(i, j, k), indexer.getIndex(i, j+1, k), 1./h));
          // tripletList.push_back(T(i + j*nx + k*nx*ny, i + (j+1)*nx + k*nx*ny, 1./h));
        }
        else{
          tripletList.push_back(T(indexer.getIndex(i, j, k), indexer.getIndex(i, j, k+1), 1./h));
          // tripletList.push_back(T(i + j*nx + k*nx*ny, i + j*nx + (k+1)*nx*ny, 1./h));
        }
      }
    }
  }

  D.setFromTriplets(tripletList.begin(), tripletList.end());
}