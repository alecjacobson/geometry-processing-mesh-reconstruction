#include "fd_partial_derivative.h"
#include <iostream>
void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  int a = 1 - ceil(dir/2.0);
  int b = ceil(dir/2.0) - floor(dir/2.0);
  int c = floor(dir/2.0);

  int new_x = nx-a;
  int new_y = ny-b;
  int new_z = nz-c;

  int next = a + (b*nx) + (c*nx*ny);

  D.resize(new_x*new_y*new_z, nx*ny*nz);

  for(int i = 0; i < new_x; i += 1){
  	for(int j = 0; j < new_y; j += 1){
  		for(int k = 0; k < new_z; k += 1){
  			D.insert(i + (new_x*j) + (new_x*new_y*k), i + (nx*j) + (nx*ny*k)) = -1;
  			D.insert(i + (new_x*j) + (new_x*new_y*k), i + (nx*j) + (nx*ny*k) + next) = 1;
  		}
  	}
  }
   std::cout << D.rows() << std::endl;
}
