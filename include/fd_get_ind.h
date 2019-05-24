#ifndef FD_GET_IND_H
#define FD_GET_IND_H
// Converts a 3D subscript to a unique index
// Input: A grid with dimensions (nx, ny, nz) and a 3D subscript (i=x, j=y, k=z)
// Output: Single integer that flattens the (i, j, k) to some index
int fd_get_ind(
  const int nx,
  const int ny,
  const int nz,
  const int x,
  const int y,
  const int z);
#endif
