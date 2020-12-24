#include "fd_partial_derivative.h"
#include <iostream>

int fast_indexing(int i, int j, int k, int nx, int ny) {
    return i + (nx * j) + (nx * ny * k);

}

void fd_partial_derivative(
        const int nx,
        const int ny,
        const int nz,
        const double h,
        const int dir,
        Eigen::SparseMatrix<double> &D) {

    int nxm = nx - (dir == 0);
    int nym = ny - (dir == 1);
    int nzm = nz - (dir == 2);


    // this triplet is <row,col,value>
    // by using this we can save time while setting the sparseM D
    typedef Eigen::Triplet<double> et;
    std::vector<et> TList;

    //number of element to insert
    TList.reserve(nxm * nym * nzm);

    D.resize(nxm * nym * nzm, nx * ny * nz);

    for (int i = 0; i < nxm; i++) {
        for (int j = 0; j < nym; j++) {
            for (int k = 0; k < nzm; k++) {
                TList.push_back(et(fast_indexing(i, j, k, nxm, nym), fast_indexing(i, j, k, nx, ny), -1));
                et n;
                switch (dir) {
                    case 0:
                        n = et(fast_indexing(i, j, k, nxm, nym), fast_indexing(i + 1, j, k, nx, ny), 1);
                        break;
                    case 1:
                        n = et(fast_indexing(i, j, k, nxm, nym), fast_indexing(i, j + 1, k, nx, ny), 1);
                        break;
                    case 2:
                        n = et(fast_indexing(i, j, k, nxm, nym), fast_indexing(i, j, k + 1, nx, ny), 1);
                        break;

                }
                TList.push_back(n);
            }
        }
    }
    D.setFromTriplets(TList.begin(), TList.end());

}
