#include "fd_interpolate.h"

class Index
{
  public:
    int nx, ny, nz;
    Index(int x, int y, int z){
        nx = x; ny = y; nz = z;
    };
    int getIndex(int i, int j, int k){
      return i + nx*j + nx*ny*k;
    };
};

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  Index indexer(nx, ny, nz);

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(P.rows() * 8);

  for (int i = 0; i < P.rows(); i++){
    double xd = (P(i, 0) - corner(0))/h;
    double yd = (P(i, 1) - corner(1))/h;
    double zd = (P(i, 2) - corner(2))/h;

    double xf = std::floor(xd);
    double yf = std::floor(yd);
    double zf = std::floor(zd);

    double xr = xd - xf;
    double yr = yd - yf;
    double zr = zd - zf;

    tripletList.push_back(T(i, indexer.getIndex(xf, yf, zf), (1-xr)*(1-yr)*(1-zr)));
    tripletList.push_back(T(i, indexer.getIndex(xf + 1, yf, zf), xr*(1-yr)*(1-zr)));
    tripletList.push_back(T(i, indexer.getIndex(xf, yf + 1, zf), (1-xr)*yr*(1-zr)));
    tripletList.push_back(T(i, indexer.getIndex(xf, yf, zf + 1), (1-xr)*(1-yr)*zr));
    tripletList.push_back(T(i, indexer.getIndex(xf + 1, yf + 1, zf), xr*yr*(1-zr)));
    tripletList.push_back(T(i, indexer.getIndex(xf + 1, yf, zf + 1), xr*(1-yr)*zr));
    tripletList.push_back(T(i, indexer.getIndex(xf, yf + 1, zf + 1), (1-xr)*yr*zr));
    tripletList.push_back(T(i, indexer.getIndex(xf + 1, yf + 1, zf + 1), xr*yr*zr));
  }

  W.setFromTriplets(tripletList.begin(), tripletList.end());
}