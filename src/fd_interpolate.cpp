#include "fd_interpolate.h"

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  double x, y, z;
  int x_floor, y_floor, z_floor;
  double x_res, y_res, z_res;
  int x_coef = 1, y_coef = 1, z_coef = 1;
  std::vector<Eigen::Triplet<double>> list;
  list.reserve(P.rows() * 8);

  for (int row = 0; row < P.rows(); row++){
    x = (P(row, 0) - corner(0)) / h;
    y = (P(row, 1) - corner(1)) / h;
    z = (P(row, 2) - corner(2)) / h;

    x_floor = std::floor(x);
    y_floor = std::floor(y);
    z_floor = std::floor(z);

    x_res = x - x_floor;
    y_res = y - y_floor;
    z_res = z - z_floor;

    for (int a = 0; a <= 1; a++){
      for (int b = 0; b <= 1; b++){
        for (int c = 0; c <= 1; c++){
          if (a == 0) x_coef = -1;
          if (b == 0) y_coef = -1;
          if (c == 0) z_coef = -1;
          list.push_back(Eigen::Triplet<double>(row,
            x_floor + a + nx * (y_floor + b) + nx * ny * (z_floor + c),
            (1 - a + x_coef * x_res) * (1 - b + y_coef * y_res) * (1 - c + z_coef * z_res)));
          x_coef = 1; y_coef = 1; z_coef = 1;
        }
      }
    }

    // list.push_back(Eigen::Triplet<double>(row, x_floor + nx * y_floor + nx * ny * z_floor, (1-x_res)*(1-y_res)*(1-z_res)));
    // list.push_back(Eigen::Triplet<double>(row, x_floor + 1 + nx * y_floor + nx * ny * z_floor, x_res*(1-y_res)*(1-z_res)));
    // list.push_back(Eigen::Triplet<double>(row, x_floor + nx * (y_floor + 1) + nx * ny * z_floor, (1-x_res)*y_res*(1-z_res)));
    // list.push_back(Eigen::Triplet<double>(row, x_floor + nx * y_floor + nx * ny * (z_floor + 1), (1-x_res)*(1-y_res)*z_res));
    // list.push_back(Eigen::Triplet<double>(row, x_floor + 1 + nx * (y_floor + 1) + nx * ny * z_floor, x_res*y_res*(1-z_res)));
    // list.push_back(Eigen::Triplet<double>(row, x_floor + 1 + nx * y_floor + nx * ny * (z_floor + 1), x_res*(1-y_res)*z_res));
    // list.push_back(Eigen::Triplet<double>(row, x_floor + nx * (y_floor + 1) + nx * ny * (z_floor + 1), (1-x_res)*y_res*z_res));
    // list.push_back(Eigen::Triplet<double>(row, x_floor + 1 + nx * (y_floor + 1) + nx * ny * (z_floor + 1), x_res*y_res*z_res));
  }
  W.resize(P.rows(), nx*ny*nz);
  W.setFromTriplets(list.begin(), list.end());
}
