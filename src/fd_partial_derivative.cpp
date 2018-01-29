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
  std::vector<Eigen::Triplet<double>> triplets;
  
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

  for(int i = 0; i < expected - off_by_one; i++){
    int neighbor = i + off_by_one;

    if (i >= off_by_one){
      // The diagonal is always 1, except when we have no neighbor 
      triplets.push_back(Eigen::Triplet<double>(i, i, 1));        
    }

    // Set neighbor ahead
    triplets.push_back(Eigen::Triplet<double>(neighbor, i, -1));        
  }

  D.setFromTriplets(triplets.begin(), triplets.end());
  cout << D.nonZeros() << endl;
}
