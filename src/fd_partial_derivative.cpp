#include "fd_partial_derivative.h"
#include <tuple>

using namespace std;

// Returns the dimensions of the staggered grid in the order (x, y, z).
tuple<int, int, int> get_staggered_grid_dimensions(
    const int nx,
    const int ny,
    const int nz,
    const int dir);

// Returns the two nodes, connected by an edge, that
// surround node i on the staggered grid. The first
// node in the tuple comes before i and the second comes after.
tuple<double, double> calculate_node_pair(
    const int i,
    const int x_dim,
    const int y_dim,
    const int nx,
    const int ny,
    const int dir);

void fd_partial_derivative(
    const int nx,
    const int ny,
    const int nz,
    const double h,
    const int dir,
    Eigen::SparseMatrix<double> & D)
{
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    tuple<int, int, int> staggered_dims = get_staggered_grid_dimensions(nx, ny, nz, dir);
    int rows = get<0>(staggered_dims) * get<1>(staggered_dims) * get<2>(staggered_dims);
    triplets.reserve(2 * rows);

    for (int i = 0; i < rows; i++)
    {
        tuple<double, double> surrounding_nodes = calculate_node_pair(i, get<0>(staggered_dims), get<1>(staggered_dims), nx, ny, dir);
        triplets.push_back(T(i, get<0>(surrounding_nodes), -h/2));
        triplets.push_back(T(i, get<1>(surrounding_nodes), h/2));
    }

    D.resize(rows, nx * ny * nz);
    D.setFromTriplets(triplets.begin(), triplets.end());
}

tuple<int, int, int> get_staggered_grid_dimensions(
    const int nx,
    const int ny,
    const int nz,
    const int dir)
{
    if (dir == X) 
    {
        return tuple<int, int, int>(nx - 1, ny, nz);
    }
    if (dir == Y) 
    {
        return tuple<int, int, int>(nx, ny - 1, nz);
    }
    
    return tuple<int, int, int>(nx, ny, nz - 1);
}

tuple<double, double> calculate_node_pair(
    const int i,
    const int x_dim,
    const int y_dim,
    const int nx,
    const int ny,
    const int dir) 
{
    int staggered_plane_size = x_dim * y_dim;
    int z_staggered = i / staggered_plane_size;
    int y_staggered = (i % staggered_plane_size) / x_dim;
    int x_staggered = (i % staggered_plane_size) % x_dim;

    double node1 = z_staggered * (nx * ny) + y_staggered * nx + x_staggered;
    double node2;

    if (dir == X) 
    {
        node2 = node1 + 1;
    }
    else if (dir == Y)
    {
        node2 = z_staggered * (nx * ny) + (y_staggered + 1) * nx + x_staggered;
    }
    else 
    {
        node2 = (z_staggered + 1) * (nx * ny) + y_staggered * nx + x_staggered;
    }

    return tuple<int, int>(node1, node2);
}
