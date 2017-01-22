#include "fd_partial_derivative.h"
#include <vector>

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
  int snx = (dir==0)?nx-1:nx;
  int sny = (dir==1)?ny-1:ny;
  int snz = (dir==2)?nz-1:nz;
  int m = snx * sny * snz;
  auto trips = fd_partial_derivative_triplets(nx,ny,nz,h,dir);
  D.resize(m,nx*ny*nz);
  D.setFromTriplets(trips.begin(),trips.end());
}
  
std::vector<Eigen::Triplet<double>> fd_partial_derivative_triplets(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir)
{
  using Triplet = Eigen::Triplet<double>;

  int snx = (dir==0)?nx-1:nx;
  int sny = (dir==1)?ny-1:ny;
  int snz = (dir==2)?nz-1:nz;
  int m = snx * sny * snz;

  auto indexer = [&](int i, int j, int k) -> int {
      return i + nx*(j + k * ny);
  };
  auto sindexer = [&](int i, int j, int k) -> int {
      return i + snx*(j + k * sny);
  };


  std::vector<Triplet> trips;
  trips.reserve(2*m);

  double C = 1.0/h;

  for(int i = 0; i < snx; ++i) {
      for(int j = 0; j < sny; ++j) {
          for(int k = 0; k < snz; ++k) {
              int si = sindexer(i,j,k);
              trips.emplace_back(Triplet{si,indexer(i,j,k),-C});
              switch(dir) {
                  case 0: trips.emplace_back(Triplet{si,indexer(i+1,j,k),C}); break;
                  case 1: trips.emplace_back(Triplet{si,indexer(i,j+1,k),C}); break;
                  case 2: trips.emplace_back(Triplet{si,indexer(i,j,k+1),C}); break;
              }
          }
      }
  }

  return trips;



}
