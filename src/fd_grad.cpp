#include "fd_grad.h"
#include "fd_partial_derivative.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
	////////////////////////////////////////////////////////////////////////////
	// Add your code here
	////////////////////////////////////////////////////////////////////////////

	Eigen::SparseMatrix<double> Gx, Gy, Gz;	
	fd_partial_derivative(nx, ny, nz, h, 0, Gx);
	fd_partial_derivative(nx, ny, nz, h, 1, Gy);
	fd_partial_derivative(nx, ny, nz, h, 2, Gz);	
	
	G.resize(Gx.rows() + Gy.rows() + Gz.rows(), nx*ny*nz);
	
	////Dosn't work this.
	//G.middleRows(0, Gx.rows()) = Gx;
	//G.middleRows(Gx.rows(), Gy.rows()) = Gy;
	//G.middleRows(Gx.rows() + Gy.rows(), Gz.rows()) = Gz;

	//Ummmm, Eigen dosn't support that way for sparse matrices.
	//Modified this example http://stackoverflow.com/questions/41756428/concatenate-sparse-matrix-eigen
	//To concatenate the gradient matrices.
	typedef Eigen::Triplet<double> Triple;
	std::vector<Triple> triplets;
	triplets.reserve(2 * (Gx.rows() + Gy.rows() + Gz.rows()));
	
	for (int i = 0; i < Gx.outerSize(); i++)
		 for (Eigen::SparseMatrix<double>::InnerIterator it(Gx, i); it; ++it)
			 triplets.push_back(Triple(it.row(), it.col(), it.value()));
	
	for (int i = 0; i < Gy.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(Gy, i); it; ++it)
			 triplets.push_back(Triple(Gx.rows() + it.row(), it.col(), it.value()));
	
	for (int i = 0; i < Gz.outerSize(); i++)
		for (Eigen::SparseMatrix<double>::InnerIterator it(Gz, i); it; ++it)
		triplets.push_back(Triple(Gx.rows() + Gy.rows() + it.row(), it.col(), it.value()));
	
	G.setFromTriplets(triplets.begin(), triplets.end());
	
	
}


