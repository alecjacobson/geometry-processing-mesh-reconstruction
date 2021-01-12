#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <iostream>

void matrix_concat_vertical(Eigen::SparseMatrix<double> & G_x, Eigen::SparseMatrix<double> & G_y,
Eigen::SparseMatrix<double> & G_z, Eigen::SparseMatrix<double> & G);

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here

  Eigen::SparseMatrix<double> G_x;
  Eigen::SparseMatrix<double> G_y;
  Eigen::SparseMatrix<double> G_z;

  fd_partial_derivative(nx, ny, nz, h, 0, G_x);
  fd_partial_derivative(nx, ny, nz, h, 1, G_y);
  fd_partial_derivative(nx, ny, nz, h, 2, G_z);
  
  matrix_concat_vertical(G_x, G_y, G_z, G);

  
  ////////////////////////////////////////////////////////////////////////////
}

void matrix_concat_vertical(Eigen::SparseMatrix<double> & G_x, Eigen::SparseMatrix<double> & G_y,
Eigen::SparseMatrix<double> & G_z, Eigen::SparseMatrix<double> & G) {

// holder for the merged sparse matrices
  std::vector<Eigen::Triplet<double>> tripletList;

  tripletList.reserve((G_x.rows() + G_y.rows() + G_z.rows())*2);

  // iterate over non-zero elements of G_x
  for (int k=0; k<G_x.outerSize(); ++k)
  for (Eigen::SparseMatrix<double>::InnerIterator it(G_x,k); it; ++it)
  {
    /*
      it.value(); // value
      it.row();   // row index
      it.col();   // col index (here it is equal to k)
      it.index(); // inner index, here it is equal to it.row()
      Ref: https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
    */

    tripletList.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
  }

  // now over G_y
  for (int k=0; k<G_y.outerSize(); ++k)
  for (Eigen::SparseMatrix<double>::InnerIterator it(G_y,k); it; ++it)
  {
    tripletList.push_back(Eigen::Triplet<double>(G_x.rows() + it.row(), it.col(), it.value()));
  }

  // now over G_z
  for (int k=0; k<G_z.outerSize(); ++k)
  for (Eigen::SparseMatrix<double>::InnerIterator it(G_z,k); it; ++it)
  {

    tripletList.push_back(Eigen::Triplet<double>(G_x.rows() + G_y.rows() + it.row(), it.col(), it.value()));
  }

  
  G.resize(G_x.rows()+G_y.rows()+G_z.rows(), G_x.cols());
  G.setZero();

  G.setFromTriplets(tripletList.begin(), tripletList.end());

}