#include "fd_interpolate.h"

#include <cmath>
#include <vector>

#include <iostream>
using std::cout;


void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
    typedef Eigen::Triplet<double> Triplets;

    Eigen::RowVector3i x;
    Eigen::RowVector3d xpos;
    std::vector<Triplets> triplets;
    triplets.reserve(P.rows() * 8);

    auto to_idx = [=](int i, int j, int k) {
        return i + nx * (j + ny * k);
    };

    for (int i = 0; i < P.rows(); ++i) {
        auto p = P.row(i);

        // Find bottom-left-front vertex `x` of the grid cell containing p
        xpos = ((p - corner) / h);
        for (int j = 0; j < 3; ++j) {
            x(j) = floor(xpos(j));
            xpos(j) = x(j) * h + corner(j);
        }

        // Weights associated with top-right-back corner of 
        //      staggered grid cell containing each point
        auto wx = (p(0) - xpos(0)) / h;
        auto wy = (p(1) - xpos(1)) / h;
        auto wz = (p(2) - xpos(2)) / h;

        triplets.emplace_back(i, to_idx(x(0), x(1), x(2)), (1-wx)*(1-wy)*(1-wz));
        triplets.emplace_back(i, to_idx(x(0)+1, x(1), x(2)), (wx)*(1-wy)*(1-wz));
        triplets.emplace_back(i, to_idx(x(0), x(1)+1, x(2)), (1-wx)*(wy)*(1-wz));
        triplets.emplace_back(i, to_idx(x(0), x(1), x(2)+1), (1-wx)*(1-wy)*(wz));
        triplets.emplace_back(i, to_idx(x(0)+1, x(1), x(2)+1), (wx)*(1-wy)*(wz));
        triplets.emplace_back(i, to_idx(x(0)+1, x(1)+1, x(2)), (wx)*(wy)*(1-wz));
        triplets.emplace_back(i, to_idx(x(0), x(1)+1, x(2)+1), (1-wx)*(wy)*(wz));
        triplets.emplace_back(i, to_idx(x(0)+1, x(1)+1, x(2)+1), (wx)*(wy)*(wz));
    }

    W.setFromTriplets(triplets.begin(), triplets.end());
}
