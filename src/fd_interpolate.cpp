#include "fd_interpolate.h"
#include <math.h>
#include "igl/floor.h"

//aux function that helps to keep the code clean and easy to follow.
int fast_indexing(Eigen::RowVector3d &p_floor, int i, int j, int k, int nx, int ny) {
    return (p_floor(0) + i) + (nx * (p_floor(1) + j)) + (nx * ny * (p_floor(2) + k));

}

void fd_interpolate(
        const int nx,
        const int ny,
        const int nz,
        const double h,
        const Eigen::RowVector3d &corner,
        const Eigen::MatrixXd &P,
        Eigen::SparseMatrix<double> &W) {


    typedef Eigen::Triplet<double> et;
    std::vector<et> TList;
    TList.reserve(8 * P.rows());


    for (int i = 0; i < P.rows(); i++) {
        // Let's normalize all the points.
        Eigen::RowVector3d p = (1. / h) * (P.row(i) - corner);
        Eigen::RowVector3d p_floor;

        // Let's compute p_0
        igl::floor(p, p_floor);

        Eigen::RowVector3d diff = p - p_floor;

        //This could be done in a for-loop but it obscures the code.
        // In any case, it is just 8 lines.
        TList.push_back(et(i, fast_indexing(p_floor, 0, 0, 0, nx, ny),
                           (1 - diff[0]) * (1 - diff[1]) * (1 - diff[2])));

        TList.push_back(et(i, fast_indexing(p_floor, 0, 0, 1, nx, ny),
                           (1 - diff[0]) * (1 - diff[1]) * (diff[2])));

        TList.push_back(et(i, fast_indexing(p_floor, 0, 1, 0, nx, ny),
                           (1 - diff[0]) * (diff[1]) * (1 - diff[2])));

        TList.push_back(et(i, fast_indexing(p_floor, 0, 1, 1, nx, ny),
                           (1 - diff[0]) * (diff[1]) * (diff[2])));

        TList.push_back(et(i, fast_indexing(p_floor, 1, 0, 0, nx, ny),
                           (diff[0]) * (1 - diff[1]) * (1 - diff[2])));

        TList.push_back(et(i, fast_indexing(p_floor, 1, 0, 1, nx, ny),
                           (diff[0]) * (1 - diff[1]) * (diff[2])));

        TList.push_back(et(i, fast_indexing(p_floor, 1, 1, 0, nx, ny),
                           (diff[0]) * (diff[1]) * (1 - diff[2])));

        TList.push_back(et(i, fast_indexing(p_floor, 1, 1, 1, nx, ny),
                           (diff[0]) * (diff[1]) * (diff[2])));
    }

    W.setFromTriplets(TList.begin(), TList.end());

}
