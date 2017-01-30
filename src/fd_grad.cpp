#include "fd_grad.h"
#include "fd_partial_derivative.h"

using namespace Eigen;

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

	// Reference http://stackoverflow.com/questions/41756428/concatenate-sparse-matrix-eigen
	std::vector<Triplet<double> > tripletList;
	for (int i = 0; i < Gx.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(Gx, i); it; ++it)
		{
			tripletList.push_back(Triplet<double>(it.row(), it.col(), it.value()));
		}
	}
	for (int i = 0; i < Gy.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(Gy, i); it; ++it)
		{
			tripletList.push_back(Triplet<double>(Gx.rows() + it.row(), it.col(), it.value()));
		}
	}
	for (int i = 0; i < Gz.outerSize(); ++i)
	{
		for (SparseMatrix<double>::InnerIterator it(Gz, i); it; ++it)
		{
			tripletList.push_back(Triplet<double>(Gx.rows() + Gy.rows()+ it.row(), it.col(), it.value()));
		}
	}

	G.resize(Gx.rows() + Gy.rows() + Gz.rows(), Gx.cols());
	G.setFromTriplets(tripletList.begin(), tripletList.end());
}