#include "poisson_surface_reconstruction.h"
#include <igl/list_to_matrix.h>
#include <igl/viewer/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>

int main(int argc, char *argv[])
{
  // Load in points + normals from .pwn file
  Eigen::MatrixXd P,N;
  {
    Eigen::MatrixXd D;
    std::vector<std::vector<double> > vD;
    std::string line;
    std::fstream in;
    in.open(argc>1?argv[1]:"../shared/data/hand.pwn");
    while(in)
    {
      std::getline(in, line);
      std::vector<double> row;
      std::stringstream stream_line(line);
      double value;
      while(stream_line >> value) row.push_back(value);
      if(!row.empty()) vD.push_back(row);
    }
    if(!igl::list_to_matrix(vD,D)) return EXIT_FAILURE;
    assert(D.cols() == 6 && "pwn file should have 6 columns");
    P = D.leftCols(3);
    N = D.rightCols(3);
  }

  // Reconstruct mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  poisson_surface_reconstruction(P,N,V,F);

  // Create a libigl Viewer object to toggle between point cloud and mesh
  igl::viewer::Viewer viewer;
  std::cout<<R"(
  P,p      view point cloud
  M,m      view mesh
)";
  const auto set_points = [&]()
  {
    viewer.data.clear();
    viewer.data.set_points(P,Eigen::RowVector3d(1,1,1));
    viewer.data.add_edges(P,(P+0.01*N).eval(),Eigen::RowVector3d(1,0,0));
  };
  set_points();
  viewer.callback_key_pressed = [&](igl::viewer::Viewer&, unsigned int key,int)
  {
    switch(key)
    {
      case 'P':
      case 'p':
        set_points();
        return true;
      case 'M':
      case 'm':
        viewer.data.clear();
        viewer.data.set_mesh(V,F);
        return true;
    }
    return false;
  };
  viewer.core.point_size = 2;
  viewer.launch();

  return EXIT_SUCCESS;
}

