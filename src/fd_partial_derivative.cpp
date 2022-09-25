#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  int x_num = nx, y_num = ny, z_num = nz;
  if (dir == 0) {
    x_num = nx - 1;
  } else if (dir == 1){
    y_num = ny - 1;
  } else {
    z_num = nz - 1;
  }

  std::vector<Eigen::Triplet<double>> list;
  for (int i = 0; i < x_num; i++){
    for (int j = 0; j < y_num; j++){
      for (int k = 0; k < z_num; k++){
        int row  = i + j * x_num + k * x_num * y_num;
        list.push_back(Eigen::Triplet<double>(row, i + nx * j + k * nx * ny, -1));

        if (dir == 0){
          list.push_back(Eigen::Triplet<double>(row, (i + 1) + j * nx + k * nx * ny, 1));
        } else if (dir == 1) {
          list.push_back(Eigen::Triplet<double>(row, i + (j + 1) * nx + k * nx * ny, 1));
        } else {
          list.push_back(Eigen::Triplet<double>(row, i + j * nx + (k + 1) * nx * ny, 1));
        }
      }
    }
  }
  D.resize(x_num * y_num * z_num, nx * ny * nz);
  D.setFromTriplets(list.begin(), list.end());
}
