#include "fd_interpolate.h"

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d &corner,
  const Eigen::MatrixXd &P,
  Eigen::SparseMatrix<double> &W) {
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here

  W.resize(P.rows(), nx * ny * nz);

  auto ind = [&](int x, int y, int z) {
    return x + nx * (y + z * ny);
  };

  // Loop through rows of P (loop through points)
  for (int i = 0; i < P.rows(); i++) {
    // Find the nearest points
    int x0 = (int) ((P(i, 0) - corner[0]) / h);
    int y0 = (int) ((P(i, 1) - corner[1]) / h);
    int z0 = (int) ((P(i, 2) - corner[2]) / h);

    // Find percent
    double px = std::fmod((P(i, 0) - corner[0]), h);
    double py = std::fmod((P(i, 1) - corner[1]), h);
    double pz = std::fmod((P(i, 2) - corner[2]), h);

    W.insert(i, ind(x0, y0, z0)) = 1 - (px / h) - (py / h) - (pz / h);
    W.insert(i, ind(x0 + 1, y0, z0)) = px / h;
    W.insert(i, ind(x0, y0 + 1, z0)) = py / h;
    W.insert(i, ind(x0, y0, z0 + 1)) = pz / h;

  }

  ////////////////////////////////////////////////////////////////////////////
}
