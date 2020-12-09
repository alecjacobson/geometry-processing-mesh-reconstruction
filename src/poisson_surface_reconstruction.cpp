#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>

#include "fd_interpolate.h"
#include "fd_partial_derivative.h"
#include "fd_grad.h"

#include <Eigen/IterativeLinearSolvers>

void poisson_surface_reconstruction(const Eigen::MatrixXd & P,
		const Eigen::MatrixXd & N, Eigen::MatrixXd & V, Eigen::MatrixXi & F) {
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
	Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad * h;
	// Grid dimensions should be at least 3
	nx = std::max(
			(P.col(0).maxCoeff() - P.col(0).minCoeff() + (2. * pad) * h) / h,
			3.);
	ny = std::max(
			(P.col(1).maxCoeff() - P.col(1).minCoeff() + (2. * pad) * h) / h,
			3.);
	nz = std::max(
			(P.col(2).maxCoeff() - P.col(2).minCoeff() + (2. * pad) * h) / h,
			3.);
	// Compute positions of grid nodes
	Eigen::MatrixXd x(nx * ny * nz, 3);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				// Convert subscript to index
				const auto ind = i + nx * (j + k * ny);
				x.row(ind) = corner + h * Eigen::RowVector3d(i, j, k);
			}
		}
	}
	Eigen::VectorXd g = Eigen::VectorXd::Zero(nx * ny * nz);

	////////////////////////////////////////////////////////////////////////////
	// Add your code here
	////////////////////////////////////////////////////////////////////////////
	// first, we calculate spaese mtx M, which M_p3_xyz3 * X_(x-1)(y-1)(z-1)3_1 = P_p3_1
	int pnum = n;
	int gnum = (nx - 1) * (ny - 1) * (nz - 1);
	Eigen::SparseMatrix<double> W(pnum * 3, gnum * 3);
	fd_interpolate(nx, ny, nz, h, corner, P, W);

	// next, we calculate sparse mtx G, which G_(x-1)(y-1)(z-1)3_xyz * X_xyz_1 = D_(x-1)(y-1)(z-1)_1
	Eigen::SparseMatrix<double> G((nx - 1) * (ny - 1) * (nz - 1) * 3,
			nx * ny * nz * 1);
	fd_grad(nx, ny, nz, h, G);

	// next, we combine G with M
	Eigen::SparseMatrix<double> A = W * G;

	// last, we solve Ax = B
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
	solver.compute(A);

	Eigen::VectorXd b(n * 3);
	for (int i = 0; i < pnum; i++) {
		for (int j = 0; j < 3; j++) {
			b(i * 3 + j, 0) = N(i, j);
		}
	}
	g = solver.solve(b);

	double sigma = (W * g).mean();
	for (int i = 0; i < g.cols(); i++)
		g(i) = g(i) - sigma;

	////////////////////////////////////////////////////////////////////////////
	// Run black box algorithm to compute mesh from implicit function: this
	// function always extracts g=0, so "pre-shift" your g values by -sigma
	////////////////////////////////////////////////////////////////////////////
	igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
