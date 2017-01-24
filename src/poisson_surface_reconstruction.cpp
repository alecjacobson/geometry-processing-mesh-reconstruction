#include "poisson_surface_reconstruction.h"
#include "fd_interpolate.h"
#include "fd_grad.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>

void poisson_surface_reconstruction(
	const Eigen::MatrixXd & P,
	const Eigen::MatrixXd & N,
	Eigen::MatrixXd & V,
	Eigen::MatrixXi & F)
{
  ////////////////////////////////////////////////////////////////////////////
  // Construct FD grid, CONGRATULATIONS! You get this for free!
  ////////////////////////////////////////////////////////////////////////////
  // number of input points
	const int n = P.rows();
	// Grid dimensions
	int nx, ny, nz;
	// Maximum extent (side length of bounding box) of points
	double max_extent =
		(P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
	  // padding: number of cells beyond bounding box of input points
	const double pad = 8;
	// choose grid spacing (h) so that shortest side gets 30+2*pad samples
	double h = max_extent / double(30 + 2 * pad);
	// Place bottom-left-front corner of grid at minimum of points minus padding
	Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad*h;
	// Grid dimensions should be at least 3 
	nx = std::max((P.col(0).maxCoeff() - P.col(0).minCoeff() + (2.*pad)*h) / h, 3.);
	ny = std::max((P.col(1).maxCoeff() - P.col(1).minCoeff() + (2.*pad)*h) / h, 3.);
	nz = std::max((P.col(2).maxCoeff() - P.col(2).minCoeff() + (2.*pad)*h) / h, 3.);
	// Compute positions of grid nodes
	Eigen::MatrixXd x(nx*ny*nz, 3);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
			   // Convert subscript to index
				const auto ind = i + nx*(j + k * ny);
				x.row(ind) = corner + h*Eigen::RowVector3d(i, j, k);
			}
		}
	}
	Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

	////////////////////////////////////////////////////////////////////////////
	// Add your code here
	////////////////////////////////////////////////////////////////////////////

	//First, construct v from W and the points P.
	//We have that P = W*x, where x is the grid locations. 
	//We want N = Wv, for each dimension.

	Eigen::VectorXd v;
	v.resize((nx - 1)*ny*nz + nx*(ny - 1)*nz + nx*ny*(nz - 1));
	{
		Eigen::SparseMatrix<double> W;
		std::vector<Eigen::Triplet<double>> v_entries;
		int row_offset = 0;
		for (int i = 0; i < 3; i++) {
			Eigen::RowVector3d offset = Eigen::RowVector3d::Zero();
			offset(i) = 1;
			fd_interpolate(nx - offset(0),
						   ny - offset(1),
						   nz - offset(2),
						   h, corner + h*offset,
						   P, W);
			auto v_temp = W.transpose()*N.col(i);
			v.block(row_offset, 0, v_temp.rows(), 1) = v_temp;
			row_offset += v_temp.rows();
		}
	}

	Eigen::SparseMatrix<double> G;
	fd_grad(nx, ny, nz, h, G);

	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
	solver.compute(G.transpose()*G);
	Eigen::VectorXd L = G.transpose()*v;
	Eigen::VectorXd solution = solver.solve(L);

	Eigen::SparseMatrix<double> W;
	fd_interpolate(nx, ny, nz, h, corner, P, W);
	double sigma = (1.0 / n) * (W * solution).sum();
	g = solution.array() - sigma;

	////////////////////////////////////////////////////////////////////////////
	// Run black box algorithm to compute mesh from implicit function: this
	// function always extracts g=0, so "pre-shift" your g values by -sigma
	////////////////////////////////////////////////////////////////////////////
	igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
