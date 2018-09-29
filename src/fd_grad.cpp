#include "fd_grad.h"
#include <vector>
#include <igl/cat.h>


void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  // Define partial derivatives for x, y, and z.
  Eigen::SparseMatrix<double> partialX, partialY, partialZ;
  fd_partial_derivative(nx, ny, nz, h, 0, partialX);
  fd_partial_derivative(nx, ny, nz, h, 1, partialY);
  fd_partial_derivative(nx, ny, nz, h, 2, partialZ);

  // From Assignment page,
  // G: (nx-1)*ny*nz+ nx*(ny-1)*nz+ nx*ny*(nz-1) by nx*ny*nz sparse gradient
  // matrix: G = [Dx;Dy;Dz]
  G.resize((nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1), nx * ny * nz);

  // Concatenate the partial derivatives to get gradient.
  // Resource: http://libigl.github.io/libigl/tutorial/
  Eigen::SparseMatrix<double> temp;
  igl::cat(1, partialX, partialY, temp);
  igl::cat(1, temp, partialZ, G);
}
