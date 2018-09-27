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

  std::vector<Eigen::Triplet<double>> tripletList;

  // Loop through rows of P (loop through points)
  for (int i = 0; i < P.rows(); i++) {
    // Find the nearest points
    int x = (int) ((P(i, 0) - corner[0]) / h);
    int y = (int) ((P(i, 1) - corner[1]) / h);
    int z = (int) ((P(i, 2) - corner[2]) / h);

    // Find percent
    double px = (P(i, 0) - corner[0]) / h - x;
    double py = (P(i, 1) - corner[1]) / h - y;
    double pz = (P(i, 2) - corner[2]) / h - z;

    tripletList.push_back(Eigen::Triplet<double> (i, ind(x, y, z), (1 - px) * (1 - py) * (1 - pz)));
    tripletList.push_back(Eigen::Triplet<double> (i, ind(x + 1, y, z), px * (1 - py) * (1 - pz)));
    tripletList.push_back(Eigen::Triplet<double> (i, ind(x, y + 1, z), (1 - px) * py * (1 - pz)));
    tripletList.push_back(Eigen::Triplet<double> (i, ind(x, y, z + 1), (1 - px) * (1 - py) * pz));
    tripletList.push_back(Eigen::Triplet<double> (i, ind(x + 1, y + 1, z), px * py * (1 - pz)));
    tripletList.push_back(Eigen::Triplet<double> (i, ind(x + 1, y, z + 1), px * (1 - py) * pz));
    tripletList.push_back(Eigen::Triplet<double> (i, ind(x, y + 1, z + 1), (1 - px) * py * pz));
    tripletList.push_back(Eigen::Triplet<double> (i, ind(x + 1, y + 1, z + 1), px * py * pz));
  }

  W.setFromTriplets(tripletList.begin(), tripletList.end());

  ////////////////////////////////////////////////////////////////////////////
}
