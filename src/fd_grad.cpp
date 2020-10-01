#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include <Eigen/Sparse>
#include <iostream>

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  ////////////////////////////////////////////////////////////////////////////
  
    std::cout << "Building gradient matrix components..." << std::flush;

	int dir;

  	// x-component
	dir = 1;
	Eigen::SparseMatrix<double> Dx;
	fd_partial_derivative(nx, ny, nz, h, dir, Dx);

  	// y-component
	dir = 2;
	Eigen::SparseMatrix<double> Dy;
	fd_partial_derivative(nx, ny, nz, h, dir, Dy);

  	// z-component
	dir = 3;
	Eigen::SparseMatrix<double> Dz;
	fd_partial_derivative(nx, ny, nz, h, dir, Dz);

	std::cout << "done." << std::endl;

	// Concatenate to create gradient
	// Note: concatenation method sourced from https://stackoverflow.com/questions/41756428/concatenate-sparse-matrix-eigen

    std::cout << "Assembling full gradient matrix..." << std::flush;

	G.resize((nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1), nx*ny*nz);
  
	// Using the triplet insertion method based on Eigen's documentation
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;

	for (int ii = 0; ii < G.cols(); ii++)
	{
	    for (Eigen::SparseMatrix<double>::InnerIterator itDx(Dx, ii); itDx; ++itDx)
	    	 tripletList.push_back(T(itDx.row(), ii, itDx.value()/h));
	         // G.insertBack(itDx.row(), ii) = itDx.value()/h;
	    
	    for (Eigen::SparseMatrix<double>::InnerIterator itDy(Dy, ii); itDy; ++itDy)
	    	 tripletList.push_back(T(itDy.row(), ii, itDy.value()/h));
	         // G.insertBack(itDy.row(), ii) = itDy.value()/h;

	    for (Eigen::SparseMatrix<double>::InnerIterator itDz(Dz, ii); itDz; ++itDz)
	    	 tripletList.push_back(T(itDz.row(), ii, itDz.value()/h));
	         // G.insertBack(itDz.row(), ii) = itDz.value()/h;
	}

    G.setFromTriplets(tripletList.begin(), tripletList.end());
	G.finalize();

	std::cout << "done." << std::endl;


  ////////////////////////////////////////////////////////////////////////////
}
