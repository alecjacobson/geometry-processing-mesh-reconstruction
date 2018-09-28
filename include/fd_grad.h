#ifndef FD_GRAD_H
#define FD_GRAD_H
#include <Eigen/Sparse>
// Construct a gradient matrix for a finite-difference grid
//
// Inputs:
//   nx  number of grid steps along the x-direction
//   ny  number of grid steps along the y-direction
//   nz  number of grid steps along the z-direction
//   h  grid step size
// Outputs:
//   G  (nx-1)*ny*nz+ nx*(ny-1)*nz+ nx*ny*(nz-1) by nx*ny*nz sparse gradient
//     matrix: G = [Dx;Dy;Dz]
//
// See also: fd_partial_derivative.h
void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G);


template <typename T>
void vstack3(
    Eigen::SparseMatrix<T>& O,
    Eigen::SparseMatrix<T>& A,
    Eigen::SparseMatrix<T>& B,
    Eigen::SparseMatrix<T>& C)
{
    typedef Eigen::Triplet<T> Triplets;
    
    std::vector<Triplets> triplets;
    triplets.reserve(A.nonZeros() + B.nonZeros() + C.nonZeros());

    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), it.value());
        }
    }

    for (int k = 0; k < B.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
            triplets.emplace_back(A.rows() + it.row(), it.col(), it.value());
        }
    }

    for (int k = 0; k < C.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
            triplets.emplace_back(A.rows() + B.rows() + it.row(), it.col(), it.value());
        }
    }

    O.setFromTriplets(triplets.begin(), triplets.end());
}


#endif
