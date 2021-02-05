#include "fd_interpolate.h"
#include <cmath>

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  double inv_h = 1/h;
  std::vector<Eigen::Triplet<double>> triples;
  triples.reserve(8 * nx * ny * nz);
  typedef Eigen::Triplet<double> T;

  for (int r = 0; r < P.rows(); r++) {
        // Find closest lattice point to P.row(r). 
        Eigen::RowVector3d current_point = P.row(r); 
        Eigen::RowVector3d dist = inv_h * (current_point - corner);
        int i = floor(dist(0));
        int j = floor(dist(1));
        int k = floor(dist(2));
        Eigen::RowVector3d lattice_point_closest = corner + (h * Eigen::RowVector3d(i, j, k));
        Eigen::RowVector3d w = inv_h * (current_point - lattice_point_closest);
        Eigen::RowVector3d inv_w = Eigen::RowVector3d::Ones() - w; 

        // Fill in interpolation values according to formula on https://en.wikipedia.org/wiki/Trilinear_interpolation
        triples.push_back(T(r, i + nx * (j + k * ny), inv_w(0) * inv_w(1) * inv_w(2))); // value at (i, j, k)
        triples.push_back(T(r, (i+1) + nx * (j + k * ny), w(0) * inv_w(1) * inv_w(2))); // value at (i+1, j, k)
        triples.push_back(T(r, i + nx * ((j+1) + k * ny), inv_w(0) * w(1) * inv_w(2))); // value at (i, j+1, k)
        triples.push_back(T(r, i + nx * (j + (k+1) * ny), inv_w(0) * inv_w(1) * w(2))); // value at (i, j, k+1)
        triples.push_back(T(r, (i+1) + nx * ((j+1) + k * ny), w(0) * w(1) * inv_w(2))); // value at (i+1, j+1, k)
        triples.push_back(T(r, (i+1) + nx * (j + (k+1) * ny), w(0) * inv_w(1) * w(2))); // value at (i+1, j, k+1)
        triples.push_back(T(r, i + nx * ((j+1) + (k+1) * ny), inv_w(0) * w(1) * w(2))); // value at (i, j+1, k+1)           
        triples.push_back(T(r, (i+1) + nx * ((j+1) + (k+1) * ny), w(0) * w(1) * w(2))); // value at (i+1, j+1, k+1)
  }
  W.setFromTriplets(triples.begin(), triples.end()); 
}
