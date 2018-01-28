#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <igl/cat.h>
#include <algorithm>
#include "fd_interpolate.h"
#include "fd_grad.h"
#include <iostream>

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
    int m = nx*ny*nz;
    int xn = m - ny*nz;
    int yn = m - nx*nz;
    int zn = m - ny*nx;

    //Construct v
    Eigen::SparseMatrix<double> Wx(n,xn),Wy(n,yn),Wz(n,zn);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    
    for(int k = 0; k<3; k++){
        std::vector<T> tripletList;
        for(int i = 0; i < n; i++){
            double coord = P(i,0) - corner(0);
            double mod = coord/h;
            int floorx = int(mod);
            double weightx = 1-(mod - floorx);
 
            coord = P(i,1) - corner(1);
            mod = coord/h;
            int floory = int(mod);
            double weighty = 1-(mod - floory);
        
            coord = P(i,2) - corner(2);
            mod = coord/h;
            int floorz = int(mod);
            double weightz = 1-(mod - floorz);
            
            if(k == 0){
                if(weightx > 0.5){
                    floorx = floorx-1;
                    weightx = weightx - 0.5;
                }
                else {
                    weightx = weightx + 0.5;
                }
            } else if(k == 1){
                if(weighty > 0.5){
                    floory = floory-1;
                    weighty = weighty - 0.5;
                }
                else {
                    weighty = weighty + 0.5;
                }
            } else {
                if(weightz > 0.5){
                    floorz = floorz-1;
                    weightz = weightz - 0.5;
                }
                else {
                    weightz = weightz + 0.5;
                }
            }
        
            tripletList.push_back(T(i,floorx + nx*floory + nx*ny*floorz,weightx*weighty*weightz));
            tripletList.push_back(T(i,floorx+1 + nx*floory + nx*ny*floorz,(1-weightx)*weighty*weightz));
            tripletList.push_back(T(i,floorx + nx*(floory+1) + nx*ny*floorz,weightx*(1-weighty)*weightz));
            tripletList.push_back(T(i,floorx+1 + nx*(floory+1) + nx*ny*floorz,(1-weightx)*(1-weighty)*weightz));
        
            tripletList.push_back(T(i,floorx + nx*floory + nx*ny*(floorz+1),weightx*weighty*(1-weightz)));
            tripletList.push_back(T(i,floorx+1 + nx*floory + nx*ny*(floorz+1),(1-weightx)*weighty*(1-weightz)));
            tripletList.push_back(T(i,floorx + nx*(floory+1) + nx*ny*(floorz+1),weightx*(1-weighty)*(1-weightz)));
            tripletList.push_back(T(i,floorx+1 + nx*(floory+1) + nx*ny*(floorz+1),(1-weightx)*(1-weighty)*(1-weightz)));
        }
        if(k == 0){
            Wx.setFromTriplets(tripletList.begin(), tripletList.end());
        } else if (k == 1) {
            Wy.setFromTriplets(tripletList.begin(), tripletList.end());
        } else {
            Wz.setFromTriplets(tripletList.begin(), tripletList.end());
        }
        tripletList.clear();
    }
    
    Eigen::VectorXd vx = (Eigen::SparseMatrix<double>(Wx.transpose()))*N.col(0);
    
    Eigen::VectorXd vy = (Eigen::SparseMatrix<double>(Wy.transpose()))*N.col(1);
    
    Eigen::VectorXd vz = (Eigen::SparseMatrix<double>(Wz.transpose()))*N.col(2);
    
    Eigen::VectorXd v(n*3,1);
    
    Eigen::VectorXd C(n*2,1);

    igl::cat(1,vx,vy,C);
    igl::cat(1,C,vz,v);
    
    //Get G
    Eigen::SparseMatrix<double> G(xn+yn+zn,m);
    fd_grad(nx,ny,nz,h,G);
    
    
//    for (int k=0; k<G.outerSize(); ++k)
//        for (Eigen::SparseMatrix<double>::InnerIterator it(G,k); it; ++it)
//        {
//            std::cout << it.value() << std::endl;
//            std::cout << it.row() << std::endl;   // row index
//            std::cout << it.col() << std::endl;   // col index (here it is equal to k)
//        }

    
    Eigen::SparseMatrix<double> Gt;
    Gt = Eigen::SparseMatrix<double>(G.transpose());
    
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.compute(Gt*G);
    g = solver.solve(Gt*v);
    
    //Get Sigma
    Eigen::SparseMatrix<double> W(n,m);
    fd_interpolate(nx,ny,nz,h,corner,P,W);
    double sigma = ((Eigen::RowVectorXd::Zero(n).array()+1).matrix()*W*g).value()/n;
    g = (g.array()-sigma).matrix();
    
//std::cout << G << std::endl;

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
