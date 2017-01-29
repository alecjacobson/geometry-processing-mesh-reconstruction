#include "fd_interpolate.h"

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
	////////////////////////////////////////////////////////////////////////////
	// Add your code here
	////////////////////////////////////////////////////////////////////////////	

	/* Comments to help me figure out wth is going on	
	Need to construct W such that P = W * x -> where P are the points and x
	are the 3D grid locations. There are (nx * ny * nz) grid positions, therefore
	the x matrix is of size (nx * ny * nz) *3. Therefore, W needs to be a matrix of 
	size (P * (nx*ny*nz)) in order for the dimensions to match. 

	Then, W represents which grid positions should be trilinearly interpolated 
	to get back the points P. The values in W are the "weights" of the interpolation
	(i.e expand the equations and isolate for each of the lattice points). Turns out
	they follow a very simple pattern.
	*/

	W.resize(P.rows(), nx * ny * nz);	// Dimension the sparse matrix

	//Use triplets instead of w.insert
	typedef Eigen::Triplet<double> Triple;
	std::vector<Triple> triplets;
	triplets.reserve(P.rows() * 8);		    

	int vx, vy, vz;
	double xd, yd, zd;

	//For each point, define a weight matrix indexed using flattened 3D indices
	for (int p = 0; p < P.rows(); p++) {
		//fprintf(stderr, "processing point %d\n", p);

		//Determine the coordinates of closest lattice block
		vx = int((P(p, 0) - corner(0)) / h);
		vy = int((P(p, 1) - corner(1)) / h);
		vz = int((P(p, 2) - corner(2)) / h);

		//Determine the points offset with relation to the lattice block
		xd = ((P(p, 0) - corner(0)) / h) - vx;
		yd = ((P(p, 1) - corner(1)) / h) - vy;
		zd = ((P(p, 2) - corner(2)) / h) - vz;

		//Populate the correct indices of the W matrix using flattened 3D indices
		triplets.push_back(Triple(p, vx + vy*nx + vz*nx*ny, (1 - xd)*(1 - yd)*(1 - zd)));	//C000
		triplets.push_back(Triple(p, vx + vy*nx + (vz + 1)*nx*ny, (1 - xd)*(1 - yd)*zd));	//C001
		triplets.push_back(Triple(p, vx + (vy + 1)*nx + vz*nx*ny, (1 - xd)*yd*(1 - zd)));	//C010
		triplets.push_back(Triple(p, vx + (vy + 1)*nx + (vz + 1) *nx*ny, (1 - xd)*yd*zd));	//C011
		triplets.push_back(Triple(p, (vx + 1) + vy*nx + vz*nx*ny, xd*(1 - yd)*(1 - zd)));	//C100
		triplets.push_back(Triple(p, (vx + 1) + vy*nx + (vz + 1)*nx*ny, xd*(1 - yd)*zd));	//C101
		triplets.push_back(Triple(p, (vx + 1) + (vy + 1)*nx + vz*nx*ny, xd*yd*(1 - zd)));	//C110
		triplets.push_back(Triple(p, (vx + 1) + (vy + 1)*nx + (vz + 1)*nx*ny, xd*yd*zd));	//C111

		//OMG my brain... okay we are done.		
	}

	W.setFromTriplets(triplets.begin(), triplets.end());

	
}
