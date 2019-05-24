#include "fd_get_ind.h"

int fd_get_ind(
  const int nx,
  const int ny,
  const int nz,
  const int x,
  const int y,
  const int z)
{
    return x + nx * (y + ny * z);
}
