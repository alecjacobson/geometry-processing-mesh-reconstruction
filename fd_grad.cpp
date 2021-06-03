#include "fd_grad.h"
#include "fd_partial_derivative.h" // used to construct Dx, Dy, and Dz

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G) {
    
    // use fd_partial_derivative to construct Dx, Dy, and Dz
    Eigen::SparseMatrix<double> Dx;
    Eigen::SparseMatrix<double> Dy;
    Eigen::SparseMatrix<double> Dz;
    
    fd_partial_derivative(nx, ny, nz, h, 0, Dx);
    fd_partial_derivative(nx, ny, nz, h, 1, Dy);
    fd_partial_derivative(nx, ny, nz, h, 2, Dz);
    
    // make G the appropriate size
    G.resize((nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1), nx*ny*nz);
    
    // build a triplet list from  the non-zero elements of Dx, Dy, and Dz
    // to be used for populating G
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(G.rows()*2);
    
    // add non-zero Dx elements to the triplet list
    for (int k=0; k<Dx.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Dx,k); it; ++it) {
            tripletList.push_back(T(it.row(), it.col(), it.value()));
        }
    }
    
    // add non-zero Dy elements to the triplet list
    for (int k=0; k<Dy.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Dy,k); it; ++it) {
            tripletList.push_back(T(it.row() + (nx-1)*ny*nz, it.col(), it.value()));
        }
    }
    
    // add non-zero Dz elements to the triplet list
    for (int k=0; k<Dz.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Dz,k); it; ++it) {
            tripletList.push_back(T(it.row() + (nx-1)*ny*nz + nx*(ny-1)*nz, it.col(), it.value()));
        }
    }
    
    // set G using the triplets list
    G.setFromTriplets(tripletList.begin(), tripletList.end());
}
