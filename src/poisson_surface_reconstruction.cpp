#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include "fd_interpolate.h"
#include "fd_grad.h"
#include<Eigen/IterativeLinearSolvers>


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
		(P.colwise().maxCoeff()-P.colwise().minCoeff()).maxCoeff();

	// padding: number of cells beyond bounding box of input points
	const double pad = 8;

	// choose grid spacing (h) so that shortest side gets 30+2*pad samples
	double h  = max_extent/double(30+2*pad);

	// Place bottom-left-front corner of grid at minimum of points minus padding
	Eigen::RowVector3d corner = P.colwise().minCoeff().array()-pad*h;

	// Grid dimensions should be at least 3 
	nx = std::max((P.col(0).maxCoeff()-P.col(0).minCoeff()+(2.*pad)*h)/h,3.);
	ny = std::max((P.col(1).maxCoeff()-P.col(1).minCoeff()+(2.*pad)*h)/h,3.);
	nz = std::max((P.col(2).maxCoeff()-P.col(2).minCoeff()+(2.*pad)*h)/h,3.);

	// Compute positions of grid nodes
	Eigen::MatrixXd x(nx*ny*nz, 3);
	for(int i = 0; i < nx; i++) 
	{
		for(int j = 0; j < ny; j++)
		{
			for(int k = 0; k < nz; k++)
			{
				// Convert subscript to index
				const auto ind = i + nx*(j + k * ny);
				x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
			}
		}
	}
	
	/// Implementation Note:
	/// Was not able to get it to correctly form the mesh. There must be a
	/// bug somewhere, but despite long hours of searching, I have not been
	/// able to solve it :(. With the exception of the bug, I am very
	/// confident in my implementation.

	// Weight matrices for each axis and the full weight matrix for the
	// non-staggered grid.
	Eigen::SparseMatrix<double> xWeightMatrix, yWeightMatrix, zWeightMatrix,
			fullWeightMatrix;
	fd_interpolate(nx - 1, ny, nz, h, corner + Eigen::RowVector3d(h / 2, 0, 0), P, xWeightMatrix);
	fd_interpolate(nx, ny - 1, nz, h, corner + Eigen::RowVector3d(0, h /  2, 0), P, yWeightMatrix);
	fd_interpolate(nx, ny, nz - 1, h, corner + Eigen::RowVector3d(0, 0, h / 2), P, zWeightMatrix);
	fd_interpolate(nx, ny, nz, h, corner, P, fullWeightMatrix);

	// Normals distributed onto grid.
	Eigen::MatrixXd vX = xWeightMatrix.transpose() * N.col(0);
	Eigen::MatrixXd vY = yWeightMatrix.transpose() * N.col(1);
	Eigen::MatrixXd vZ = zWeightMatrix.transpose() * N.col(2);
	Eigen::MatrixXd v(vX.rows() + vY.rows() + vZ.rows(), 1);
	v << vX, vY, vZ;

	// Gradient matrix via finite differences.
	Eigen::SparseMatrix<double> gradient;
	fd_grad(nx, ny, nz, h, gradient);

	// Solve for the grid values. The standard form for least squares is 
	// something like: Ax = b, so we have to compute the square gradient and
	// apply the gradient transpose to v.
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	Eigen::SparseMatrix<double> squareGradient = gradient.transpose() * gradient;
	solver.compute(squareGradient);
	Eigen::VectorXd g = solver.solve(gradient.transpose() * v);

	// Find the iso value and shift g by it.
	Eigen::MatrixXd sigmaProduct = fullWeightMatrix * g;
	printf("Sigma product size: %d, %d\n", sigmaProduct.rows(), sigmaProduct.cols());
	double sigma = sigmaProduct.sum() / n;
	printf("Sigma: %f\n", sigma);
	g.array() -= sigma;

	////////////////////////////////////////////////////////////////////////////
	// Run black box algorithm to compute mesh from implicit function: this
	// function always extracts g=0, so "pre-shift" your g values by -sigma
	////////////////////////////////////////////////////////////////////////////
	igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);

	printf("Done\n");
}
