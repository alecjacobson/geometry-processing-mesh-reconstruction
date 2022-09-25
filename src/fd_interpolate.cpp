#include "fd_interpolate.h"
#include <tuple>

double get_col_index_from_coord(Eigen::Vector3d v, double h, const Eigen::RowVector3d & corner, int nx, int ny)
{
    Eigen::RowVector3d w = (v.transpose() - corner) * 1/h;
    Eigen::RowVector3d u(round(w(0)), round(w(1)), round(w(2)));
    return u(2) * (nx * ny) + u(1) * nx + u(0);
}

std::tuple<double, double> linearly_interpolate(double x0, double x1, double p) 
{
    double w0 = (p - x1) / (x0 - x1);
    double w1 = 1 - w0;
    return std::tuple<double, double>(w0, w1);
}

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(8 * P.rows());

    for (int i = 0; i < P.rows(); i++) 
    {
        Eigen::Vector3d p = P.row(i);

        Eigen::Vector3d x0(floor((p(0) - corner(0)) / h) * h + corner(0), floor((p(1) - corner(1)) / h) * h + corner(1), floor((p(2) - corner(2)) / h) * h + corner(2));
        Eigen::Vector3d x1(x0(0) + h, x0(1), x0(2));
        Eigen::Vector3d x2(x0(0), x0(1) + h, x0(2));
        Eigen::Vector3d x3(x0(0), x0(1), x0(2) + h);
        Eigen::Vector3d x4(x0(0) + h, x0(1) + h, x0(2));
        Eigen::Vector3d x5(x0(0) + h, x0(1), x0(2) + h);
        Eigen::Vector3d x6(x0(0), x0(1) + h, x0(2) + h);
        Eigen::Vector3d x7(x0(0) + h, x0(1) + h, x0(2) + h);

        std::tuple<double, double> t0 = linearly_interpolate(x0(0), x1(0), p(0));
        double w0 = std::get<0>(t0);
        double w1 = std::get<1>(t0);

        Eigen::Vector3d y0(w0 * x0(0) + w1 * x1(0), x0(1), x0(2));
        Eigen::Vector3d y1(w0 * x2(0) + w1 * x4(0), x2(1), x2(2));
        Eigen::Vector3d y2(w0 * x3(0) + w1 * x5(0), x3(1), x3(2)); 
        Eigen::Vector3d y3(w0 * x6(0) + w1 * x7(0), x6(1), x6(2));

        std::tuple<double, double> t1 = linearly_interpolate(y0(1), y1(1), p(1));
        double w2 = std::get<0>(t1);
        double w3 = std::get<1>(t1);
        
        Eigen::Vector3d z0(y0(0), w2 * y0(1) + w3 * y1(1), y0(2));
        Eigen::Vector3d z1(y2(0), w2 * y2(1) + w3 * y3(1), y2(2));

        std::tuple<double, double> t2 = linearly_interpolate(z0(2), z1(2), p(2));
        double w4 = std::get<0>(t2);
        double w5 = std::get<1>(t2);

        triplets.push_back(T(i, get_col_index_from_coord(x0, h, corner, nx, ny), w4 * w2 * w0));
        triplets.push_back(T(i, get_col_index_from_coord(x1, h, corner, nx, ny), w4 * w2 * w1));
        triplets.push_back(T(i, get_col_index_from_coord(x2, h, corner, nx, ny), w4 * w3 * w0));
        triplets.push_back(T(i, get_col_index_from_coord(x3, h, corner, nx, ny), w5 * w2 * w0));
        triplets.push_back(T(i, get_col_index_from_coord(x4, h, corner, nx, ny), w4 * w3 * w1));
        triplets.push_back(T(i, get_col_index_from_coord(x5, h, corner, nx, ny), w5 * w2 * w1));
        triplets.push_back(T(i, get_col_index_from_coord(x6, h, corner, nx, ny), w5 * w3 * w0));
        triplets.push_back(T(i, get_col_index_from_coord(x7, h, corner, nx, ny), w5 * w3 * w1));
    }

    W.resize(P.rows(), nx * ny * nz);
    W.setFromTriplets(triplets.begin(), triplets.end());
}