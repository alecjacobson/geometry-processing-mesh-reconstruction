#include "fd_partial_derivative.h"
#include <iostream>

using namespace std;

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  cout << "expecting non-zeros: " << 2 * D.rows() << endl;
  int num_rows_D = D.rows();

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(2 * num_rows_D);
  
  int expected;
  if (dir == 0){
    expected = (nx - 1) * ny * nz;
  }
  if (dir == 1){
    expected = nx * (ny - 1) * nz;
  }
  if (dir == 2){
    expected = nx * ny * (nz - 1);
  }

  int off_by_one;
  if (dir == 0){
    off_by_one = 1;
  }
  if (dir == 1){
    off_by_one = nx;
  }
  if (dir == 2){
    off_by_one = nx * ny;
  }  

  for(int i = off_by_one; i < expected + off_by_one; i++){
    // We are reasoning about the physical point i, but the matrix is missing the first example

    // Diagonal is always -1
    triplets.push_back(Eigen::Triplet<double>(i - off_by_one, i - off_by_one, -1));     
    
    // In the space we are interested,  
    triplets.push_back(Eigen::Triplet<double>(i - off_by_one, i, 1));        
  }

  D.setFromTriplets(triplets.begin(), triplets.end());
  cout << "got non-zeros: " << D.nonZeros() << endl;
}
