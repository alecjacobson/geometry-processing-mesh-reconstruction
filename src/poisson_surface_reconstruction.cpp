#include "poisson_surface_reconstruction.h"
#include "fd_interpolate.h"
#include "fd_partial_derivative.h"
#include "fd_grad.h"

#include <Eigen/SparseCore>
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <iostream> // for debugging

void testInterpolation();

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
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

//std::cout << "nx,ny,nz " << nx << ", " << ny << ", " << nz << std::endl;

  //std::cout << "testInterpolation..." << std::endl;
  //testInterpolation();
  
  // std::cout << "Setting up partial_derivative call for debugging only..." << std::endl;
  // Eigen::SparseMatrix<double> Dx, Dy, Dz;
  // fd_partial_derivative( nx, ny, nz, h, 0, Dx );
  // fd_partial_derivative( nx, ny, nz, h, 1, Dy );
  // fd_partial_derivative( nx, ny, nz, h, 2, Dz );

  // std::cout << "Setting up partial_derivative call for debugging only..." << std::endl;
  // Eigen::SparseMatrix<double> Dx, Dy, Dz;
  // fd_partial_derivative( 3, 3, 3, 1, 0, Dx );
  // fd_partial_derivative( 3, 3, 3, 1, 1, Dy );
  // fd_partial_derivative( 3, 3, 3, 1, 2, Dz );
  // return;
  
  // std::cout << "Setting up grad call for debugging only..." << std::endl;
  // Eigen::SparseMatrix<double> testG;
  // //fd_grad( nx, ny, nz, h, G );
  // fd_grad( 3, 3, 3, 1, testG );
  // std::cout<<"testG is:"<<std::endl<<testG<<std::endl;
  // return;
                                                               

  Eigen::SparseMatrix<double> W;
  fd_interpolate( nx, ny, nz, h, corner, P, W );
  
  // remember, staggered grids are offset 1/2 from the primary!
  // that means the corner too!
  Eigen::RowVector3d X(1,0,0), Y(0,1,0), Z(0,0,1);
  const double hh = 0.5 * h;
  
  Eigen::SparseMatrix<double> Wx, Wy, Wz;
  fd_interpolate( nx-1, ny, nz, h, corner + ( hh * X ), P, Wx );
  fd_interpolate( nx, ny-1, nz, h, corner + ( hh * Y ), P, Wy );
  fd_interpolate( nx, ny, nz-1, h, corner + ( hh * Z ), P, Wz );
  
  // distribute normals onto grids for V
  //V.resize( Wx.cols() + Wy.cols() + Wz.cols() );
  // x
  Eigen::VectorXd Vx = Wx.transpose() * N.col(0);
  Eigen::VectorXd Vy = Wy.transpose() * N.col(1);
  Eigen::VectorXd Vz = Wz.transpose() * N.col(2);

  std::cout << "Size of Vx is " << Vx.rows() << "x" << Vx.cols() << std::endl;
  std::cout << "Size of Vy is " << Vy.rows() << "x" << Vy.cols() << std::endl;
  std::cout << "Size of Vz is " << Vz.rows() << "x" << Vz.cols() << std::endl;
  std::cout << "Size of Wx is " << Wx.rows() << "x" << Wx.cols() << std::endl;
  std::cout << "Size of Wy is " << Wy.rows() << "x" << Wy.cols() << std::endl;
  std::cout << "Size of Wz is " << Wz.rows() << "x" << Wz.cols() << std::endl;
  //std::cout << Vx;
  //V.resize( Wx.cols() + Wy.cols() + Wz.cols() );
  Eigen::VectorXd v( Wx.cols() + Wy.cols() + Wz.cols() );
  //stack V components to make V
  v.segment( 0, Wx.cols() )                     = Vx;
  v.segment( Wx.cols(), Wy.cols() )             = Vy;
  v.segment( Wx.cols() + Wy.cols(), Wz.cols() ) = Vz;

  std::cout << "Size of v is " << v.rows() << "x" << v.cols() << std::endl;
  std::cout << "Check sum(v): " << v.sum() << std::endl;

  // Now build G from Dx, Dy, Dz by stacking
  Eigen::SparseMatrix<double> G;
  fd_grad( nx, ny, nz, h, G );
  // std::vector< Eigen::Triplet<double> > t;
  // extractTriplets( Dx, t );
  // extractTriplets( Dy, t );
  // extractTriplets( Dz, t );
  // G.setFromTriplets( t.begin(), t.end() );
  // G.finalize();

  // DEBUGGING HERE
  //std::cout << "ans:" << G.transpose() * Eigen::VectorXd::Ones(G.rows()) << std::endl;
  //return;
  
  // G^T G g = G^T v
  // rewrite as -- G2 g = Gv (Ax = b)
  Eigen::SparseMatrix<double> G2 = G.transpose() * G;
  Eigen::VectorXd Gv = G.transpose() * v;

  std::cout << "Solving G^T G g = G^T v in the form of Ax = b for unknown g" << std::endl;
  std::cout << "Solving with size of A " << G2.rows() << "x" << G2.cols() << std::endl;
  std::cout << "Solving with size of B " << Gv.rows() << "x" << Gv.cols() << std::endl;
  

  // Try different solvers...
  //Eigen::BiCGSTAB< Eigen::SparseMatrix<double> 
  Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > cgSolver;
  std::cout << "Computing G^T G..." << std::endl;
  cgSolver.compute( G2 );
  std::cout << "done..." << std::endl;
  std::cout << "Solving for g..." << std::endl;
  g = cgSolver.solve( Gv );
  std::cout << "done..." << std::endl;

  std::cout << "Computing sigma...the average of the values." << std::endl;
  double sigma = ( W * g ).array().sum() / W.rows();
  g.array() -= sigma;
  std::cout << "sigma was: " << sigma << std::endl;
  
  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}


void
testInterpolation()
{
  std::cout << "Setting up fd_interpolate call for debuging only..." << std::endl;
  Eigen::SparseMatrix<double> W;

  // Let's have 1 P
  Eigen::MatrixXd P(1,3);
  P << .25, .25, .25;
  //      .5,  .5,   .5;
  // P << 0,0,0,
  //      1,1,1;
  Eigen::RowVector3d corner = {0,0,0}; //origin
  fd_interpolate( 3, 3, 3, 1, corner, P, W );
  std::cout << "Points P= " << std::endl << P << std::endl;
  std::cout << "W is a " << W.rows() << "x" << W.cols() << std::endl;
  std::cout << "W= " << std::endl << W << std::endl;
  return;
    
}
