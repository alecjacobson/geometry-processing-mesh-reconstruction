#include "poisson_surface_reconstruction.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include "fd_grad.h"
#include "fd_interpolate.h"

void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F) {
    
    // Construct FD grid
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
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            for(int k = 0; k < nz; k++){
                // Convert subscript to index
                const auto ind = i + nx*(j + k * ny);
                x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
            }
        }
    }
    
    Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);
    
    // build gradient matrix for finite difference grid
    Eigen::SparseMatrix<double> G;
    fd_grad(nx, ny, nz, h, G);

    // build trilinear interpolation weight matrix for the regular grid
    Eigen::SparseMatrix<double> W;
    fd_interpolate(nx, ny, nz, h, corner, P, W);
   
    // build the trilinear interpolation weight matrix for the staggered grids
    Eigen::SparseMatrix<double> Wx;
    Eigen::SparseMatrix<double> Wy;
    Eigen::SparseMatrix<double> Wz;
    
    // to interpolate over the staggered grids, shift the corner by h/2 in the appropriate
    // axis and decrement nx, ny, or nz
    fd_interpolate(nx-1, ny, nz, h, corner + Eigen::RowVector3d(h/2,0,0), P, Wx);
    fd_interpolate(nx, ny-1, nz, h, corner + Eigen::RowVector3d(0,h/2,0), P, Wy);
    fd_interpolate(nx, ny, nz-1, h, corner + Eigen::RowVector3d(0,0,h/2), P, Wz);
    
    // distribute the normal values onto the corresponding staggered grid and concatenate them
    // to assemble v 
    Eigen::VectorXd vx = Wx.transpose()*(N.col(0));
    Eigen::VectorXd vy = Wy.transpose()*(N.col(1));
    Eigen::VectorXd vz = Wz.transpose()*(N.col(2));

    Eigen::VectorXd v(vx.size() + vy.size() + vz.size());
    v << vx, vy, vz;
    
    // build and solve the linear system 
    Eigen::VectorXd b;
    Eigen::SparseMatrix<double> A;
    
    A = G.transpose()*G;
    b = G.transpose()*v;

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    
    solver.compute(A);
    g = solver.solve(b);
    
    // determine a good iso-level to extract from g
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd sigma = Eigen::VectorXd::Constant(g.size(), (1.0/n)*(ones.dot(W*g)));

    // pre-shift g values by -sigma
    g -= sigma;
    
    // Black box
    igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
