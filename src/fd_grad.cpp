#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include "igl/cat.h"

void fd_grad(
        const int nx,
        const int ny,
        const int nz,
        const double h,
        Eigen::SparseMatrix<double> &G) {


    Eigen::SparseMatrix<double> Dx, Dy, Dz, buff;
    fd_partial_derivative(nx, ny, nz, h, 0, Dx);
    fd_partial_derivative(nx, ny, nz, h, 1, Dy);
    fd_partial_derivative(nx, ny, nz, h, 2, Dz);

    //matrix: G = [Dx;Dy;Dz]

    igl::cat(1, Dx, Dy, buff);
    igl::cat(1, buff, Dz, G);

}
