#include "fd_interpolate.h"

typedef Eigen::Triplet<double> tri;

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  // reserve the amount of space needed for the sparse matrix (8 nearest grid corners to the given point)
  std::vector<tri> values;
  values.reserve(8*P.rows());
  W.resize(P.rows(), nx*ny*nz);
  
  Eigen::RowVector3d gridLocation;
  int32_t cornerX = 0, cornerY = 0, cornerZ = 0;
  double  localX  = 0, localY  = 0, localZ  = 0;
    
  for( int32_t i = 0 ; i < P.rows() ; i++ )
  {
    // get the grid coords of the point:
    gridLocation = (P.row(i) - corner)/h;

    // trilinear interpolation in this context:  Wv = n
    // scale the space into unit cube to make calcs easier and then form interp weights
    // front corner locations
    cornerX = floor(gridLocation.x());
    cornerY = floor(gridLocation.y());
    cornerZ = floor(gridLocation.z());

    // local coords (within cell)
    localX = gridLocation.x() - cornerX;
    localY = gridLocation.y() - cornerY;
    localZ = gridLocation.z() - cornerZ;

    // push the values into the W matrix:
    for( int32_t x = 0 ; x < 2 ; x++ )
    {
      for( int32_t y = 0 ; y < 2 ; y++ )
      {
	for( int32_t z = 0 ; z < 2 ; z++ )
	{
	  // using some tricky indexing here
	  values.push_back(tri(i, cornerX+x + nx*(cornerY + y) + nx*ny*(cornerZ + z),
			       (x*(localX) + (1-x)*(1-localX))*(y*(localY) + (1-y)*(1-localY))*(z*(localZ) + (1-z)*(1-localZ)) ));
	}
      }
    }
    
  }

  W.setFromTriplets(values.begin(), values.end());
  
}
