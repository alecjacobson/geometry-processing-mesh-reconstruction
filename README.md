# Geometry Processing – Mesh Reconstruction

> **To get started:** Fork this repository then issue
> 
>     git clone --recursive http://github.com/[username]/geometry-processing-mesh-reconstruction.git
>

## Installation, Layout, and Compilation

See
[introduction](http://github.com/[username]/geometry-processing-introduction).

## Execution

Once built, you can execute the assignment from inside the `build/` using 

    ./mesh-reconstruction [path to point cloud]

## Background

In this assignment, we will be implementing a simplified version of the method
in  "Poisson Surface Reconstruction" by Kazhdan et al. 2006. (Your first "task"
will be to read and understand this paper).

Many scanning technologies output a set of $n$ point samples $\P$ on the 
surface of the object in question. From these points and perhaps the location
of the camera, one can also estimate normals $\N$ to the surface for each point
$\p ∈ \P$.

For shape analysis, visualization and other downstream geometry processing
phases, we would like to convert this finitely sampled _point cloud_ data into
an _explicit continuous surface representation_: i.e., a [triangle
mesh](https://en.wikipedia.org/wiki/Triangulation%5F(topology)) (a special case
of a [polygon mesh](https://en.wikipedia.org/wiki/Polygon%5Fmesh)).

### Voxel-based Implicit Surface

Converting the point cloud directly to a triangle mesh makes it very difficult
to ensure that the mesh meets certain _topological_ postconditions: i.e., that
it is [manifold](https://en.wikipedia.org/wiki/Piecewise%5Flinear%5Fmanifold),
[closed](https://en.wikipedia.org/wiki/Manifold#Manifold%5Fwith%5Fboundary),
and has a small number of
[holes](https://en.wikipedia.org/wiki/Genus%5F(mathematics)).

Instead we will first convert the point cloud _sampling representation_ into a
an _[implicit surface
representation](https://en.wikipedia.org/wiki/Implicit%5Fsurface)_: where the
unknown surface is defined as the
[level-set](https://en.wikipedia.org/wiki/Level%5Fset) of some function $g: \R³
→ \R$ mapping all points in space to a scalar value. For example, we may define
the surface $∂\S$ of some solid, volumetric shape $\S$ to be all points $\x ∈
\R³$ such that $g(x) = σ$, where we may arbitrarily set $σ=½$.

\\[
∂\S = \{\x ∈ \R³ | g(\x)  = σ\}.
\\]

On the computer, it is straightforward
[discretize](https://en.wikipedia.org/wiki/Discretization) an implicit
function. We define a regular 3D grid of
[voxels](http://en.wikipedia.org/wiki/Voxel) containing at least the [bounding
box](https://en.wikipedia.org/wiki/Minimum%5Fbounding%5Fbox) of $\S$. At each
node in the grid $\x_{i,j,k}$ we store the value of the implicit function
$g(\x_{i,j,k})$. This defines $g$ _everywhere_ in the grid via [trilinear
interpolation](https://en.wikipedia.org/wiki/Trilinear_interpolation). 

![](images/trilinear-interpolation.jpg)

For example, consider a point $\x = (⅓,¼,½)$ lying in the middle of the
bottom-most, front-most, left-most cell. We know the values at the eight
corners. Trilinear interpolation can be understood as [linear
interpolation](https://en.wikipedia.org/wiki/Linear_interpolation) in the
$x$-direction by $⅓$ on each $x$-axis-aligned edge, resulting in four values
_living_ on the same plane. These can then be linearly interpolated in the $y$
direction by $¼$ resulting in two points on the same line, and finally in the
$z$ direction by $½$ to get to our evaluation point $(⅓,¼,½)$.

An implicit surface stored as the level-set of a trilinearly interpolated grid
can be _contoured_ into a triangle mesh via the [Marching Cubes
Algorithm](https://en.wikipedia.org/wiki/Marching&5Fcubes). 
For the purposes of this assignment, we will treat this as a [black
box](https://en.wikipedia.org/wiki/Black%5Fbox). Instead, we focus on
determining what values for $g$ to store on the grid.

## Characteristic functions of solids

We assume that our set of points $\P$ lie on the surface $∂\S$ of some physical
[solid](https://en.wikipedia.org/wiki/Solid) object $\S$. This solid object
must have some non-trivial volume that we can calculate _abstractly_ as the
integral of unit density over the solid:

\\[
∫\limits_\S 1 \;dA.                                                        %_
\\]

We can rewrite this definite integral as an indefinite integral over all of
$\R^3$:

\\[
∫\limits_{\R³} χ_\S(\x) \;dA,
\\]

by introducing the [characteristic
function](https://en.wikipedia.org/wiki/Indicator%5Ffunction) of $\S$, that is
_one_ for points inside of the shape and _zero_ for points outside of $\S$:

\\[
χ_\S(\x) = \begin{cases}
  1 & \text{ if $\x ∈ \S$ } \\
  0 & \text{ otherwise}.
\end{cases}                                                             %_
\\]

Compared to typical [implicit surface
functions](https://en.wikipedia.org/wiki/Implicit%5Fsurface), this function
represents the surface $∂\S$ of the shape $\S$ as the _discontinuity_ between
the one values and the zero values. Awkwardly, the gradient of the
characteristic function $∇χ_\S$ is _not defined_ along $∂\S$.

One of the key observations made in [Kazhdan et al. 2006] is that the gradient
of a infinitesimally [mollified](https://en.wikipedia.org/wiki/Mollifier)
(smoothed) characteristic function: 

  1. points in the direction of the normal near the surface $∂\S$, and 
  2. is zero everywhere else.

Our goal will be to use our points $\P$ and normals $\N$ to _optimize_ an
implicit function $g$ over a regular grid, so that its gradient $∇g$ meets
these two criteria. In that way, our $g$ will be an approximation of the
mollified characteristic function.

## Poisson surface reconstruction

### Or: how I learned to stop worrying and minimize squared Gradients

<span id=assumptions>
Let us start by making two assumptions:
</span>

  1. we know how to compute $∇g$ at each node location $\x_{i,j,k}$, and
  2. our input points $\P$ all lie perfectly at grid nodes: 
  $∃\ \x_{i,j,k} = \p_ℓ$.

(We will find out these assumptions are not realistic, but let's defer that
discussion for now.)

If our points $\P$ lie at grid points, then our corresponding normals $\N$ also
_live_ at grid points. This leads to a very simple set of linear equations to
define a function $g$ with a gradient equal to the surface normal at the
surface and zero gradient away from the surface:

\\[
∇g(\x_{i,j,k}) = \v_{i,j,k} := \begin{cases}
  \vphantom{\left(\begin{array}{c}0\\0\\0\end{array}\right)}
  \n_ℓ & \text{ if $∃\ \p_ℓ = \x_{i,j,k}$}, \\
  \left(\begin{array}{c}
    0\\
    0\\
    0\end{array}\right) & \text{ otherwise}.
\end{cases}
\\]

This is a _vector-valued_ equation. The gradients, normals and zero-vectors are
three-dimensional (e.g., $∇g ∈ \R³$). In effect, this is _three equations_ for
every grid node.

Since we only need a single number at each grid node (the value of $g$), we
have _too many_ equations.

Like many geometry processing algorithms confronted with such an [over
determined](https://en.wikipedia.org/wiki/Overdetermined%5Fsystem), we will
_optimize_ for the solution that best _minimizes_ the error of equation:

\\[
‖ ∇g(\x_{i,j,k})  - \v_{i,j,k}‖².
\\]

We will treat the error of each grid location equally by minimizing the sum
over all grid locations:

\\[
\min_\g ∑_i ∑_j ∑_k ½ ‖ ∇g(\x_{i,j,k})  - \v_{i,j,k}‖²,
\\]

where $\g$ (written in boldface) is a vector of _unknown_ grid-nodes values,
where $g_{i,j,k} = g(\x_{i,j,k})$. 

Part of the convenience of working on a regular grid is that we can use the
[finite difference
method](https://en.wikipedia.org/wiki/Finite_difference_method) to approximate
the gradient $∇g$ on the grid.

After revisiting [our assumptions](#assumptions), we will be able to compute
approximations of 
the $x$-, $y$- and $z$-components of $∇g$ via a [sparse
matrix](https://en.wikipedia.org/wiki/Sparse%5Fmatrix) multiplication of
a "gradient matrix" $\G$ and our vector of unknown grid values $\g$. We will be
able to write the
minimization problem above in matrix form:

\\[
\min_\g ½ ‖ \G \g - \v ‖²,
\\]

or equivalently after expanding the norm:

\\[
\min_\g ½ \g^\transpose \G^\transpose \G \g - \g^\transpose \G^\transpose \v + \text{constant},
\\]

This is a quadratic "energy" function of the variables of $\g$, its minimum occurs when
an infinitesimal change in $\g$ produces no change in the energy:

\\[
\frac{∂}{∂\g} ½ \g^\transpose \G^\transpose \G \g - \g^\transpose \G^\transpose \v = 0.
\\]

Applying this derivative gives us a _sparse_ system of linear equations

\\[
\G^\transpose \G \g = \G^\transpose \v.
\\]

We will assume that we can solve this using a black box sparse solver.

Now, let's revisit [our assumptions](#assumptions).

### Gradients on a regular grid

The gradient of a function $g$ in 3D is nothing more than a vector containing
partial derivatives in each coordinate direction:

\\[
∇g(\x) = \left(\begin{array}{c}
    \frac{∂g(\x)}{∂x} \\
    \frac{∂g(\x)}{∂y} \\
    \frac{∂g(\x)}{∂z}
  \end{array}\right).
\\]

We will approximate each partial derivative individually. Let's consider the
partial derivative in the $x$ direction, $∂g(\x)/∂x$, and we will assume
without loss of generality that what we derive applies _symmetrically_ for $y$
and $z$.

The partial derivative in the $x$-direction is a one-dimensional derivative.
This couldn't be easier to do with finite differences. We approximate the
derivative of the function $g$ with respect to the $x$ direction is the
difference between the function evaluated at one grid node and at the grid node
_before_ it in the $x$-direction then divided by the spatial distance between
adjacent nodes $h$ (i.e., the grid step size):

\\[
\frac{∂g(\x_{i-½,j,k})}{∂x} = \frac{g_{i,j,k} - g_{i-1,j,k}}{h},
\\]

where we use the index $i-½$ to indicate that this derivative in the
$x$-direction lives on a [staggered
grid](https://en.wikipedia.org/wiki/Staggered%5Fgrid) _in between_ the grid
nodes where the function values for $g$.

The following pictures show a 2D example, where $g$ lives on the nodes of a
5×5 blue grid:

![](images/primary-grid.jpg)

The partial derivatives of $g$ with respect to the $x$-direction $∂g(\x)/∂x$
live on a **4**×5 green, staggered grid:

![](images/staggered-grid-x.jpg)

The partial derivatives of $g$ with respect to the $y$-direction $∂g(\x)/∂y$
live on a 5×**4** yellow, staggered grid:

![](images/staggered-grid-x-and-y.jpg)

Letting $\g ∈ \R^{n_xn_yn_z × 1}$ be column vector of function values on the
_primary grid_ (blue in the example pictures), we can construct a sparse matrix
$\D^x ∈ \R^{(n_x-1)n_yn_z × n_xn_yn_z}$ so that each row $\D^x_{i-½,j,k} ∈
\R^{1 × n_xn_yn_z}$ computes the partial derivative at the corresponding
staggered grid location $\x_{i-½,j,k}$. The $ℓ$th entry in that row receives a
value only for neighboring primary grid nodes:

\\[
\D^x_{i-½,j,k}(ℓ) = 
  \begin{cases}
  -1 & \text{ if $ℓ = i-1$ }\\
   1 & \text{ if $ℓ = i$ }\\
   0 & \text{ otherwise}
  \end{cases}.
\\]

> #### Indexing 3D arrays
> 
> Now, obviously in our code we cannot _index_ the column vector $\g$ by a
> triplet of numbers $\{i,j,k\}$ or the rows of $\D^x$ by the triplet
> ${i-½,j,k}$. We will assume that $\g_{i,j,k}$ refers to
> `g(i+j*n_x+k*n_y*n_x)`. Similarly, for the staggered grid subscripts
> ${i-½,j,k}$ we will assume that $\D^x_{i-½,j,k}(ℓ)$ refers to the matrix
> entry `Dx(i+j*n_x+k*n_y*n_x,l)`, where the $i-½$ has been _rounded down_.
>

We can similarly build matrices $\D^y$ and $\D^z$ and _stack_ these matrices
vertically to create a gradient matrix $\G$:

\\[
\G = 
\left(\begin{array}{c}
  \D^x \\
  \D^y \\
  \D^z
\end{array}\right)
  ∈ \R^{ \left((n_x-1)n_yn_z + n_x(n_y-1)n_z + n_xn_y(n_z-1)\right)×n_xn_yn_z}
\\]

This implies that our vector $\v$ of zeros and normals in our minimization
problem should not _live_ on the primary, but rather it, too, should be broken
into $x$-, $y$- and $z$-components that live of their resepctive staggered
grids:

\\[
\v = 
\left(\begin{array}{c}
  \v^x \\
  \v^y \\
  \v^z 
\end{array}\right)
  ∈ \R^{ \left((n_x-1)n_yn_z + n_x(n_y-1)n_z + n_xn_y(n_z-1)\right)×1}.
\\]

This leads to addressing our second assumption.

### B-b-b-b-but the input normals might not be at grid node locations?

At this point, we would _actually_ liked to have had that our input normals
were given component-wise on the staggered grid. Then we could immediate stick
them into $\v$. But this doesn't make much sense as each normal $\n_ℓ$ _lives_
at its associated point $\p_ℓ$, regardless of any grids.

To remedy this, we will distribute each component of each input normal $\n_ℓ$
to $\v$ at the corresponding staggered grid node location.

For example, consider the normal $\n$ at some point $\x_{1,¼,½}$. Conceptually,
we'll think of the $x$-component of the normal $n_x$ as floating in the
staggered grid corresponding to $\D^x$, in between the eight staggered grid
locations:

\\[
\x_{½,0,0},  \\
\x_{1½,0,0},   \\
\x_{½,1,0},  \\
\x_{1½,1,0},   \\
\x_{½,0,1},  \\
\x_{1½,0,1},   \\
\x_{½,1,1}, \text{ and } \\
\x_{1½,1,1}
\\]

Each of these staggered grid nodes has a corresponding $x$ value in the vector
$\v^x$.

We will distribute $n_x$ to these entries in $\v^x$ by _adding_ a partial amount
of $n_x$ to each. I.e., 

\\[
v^x_{½,0,0}  = w_{ ½,0,0}\left(\x_{1,¼,½}\right)\ n_x,  \\
v^x_{1½,0,0} = w_{1½,0,0}\left(\x_{1,¼,½}\right)\ n_x,  \\
\vdots \\
v^x_{1½,1,1} = w_{1½,1,1}\left(\x_{1,¼,½}\right)\ n_x.
\\]
where $w_{iِ+½,j,k}(\p)$ is the trilinear interpolation _weight_ associate with
staggered grid node $\x_{iِ+½,j,k}$ to interpolate a value at the point $\p$.
The trilinear interpolation weights so that:

\\[
n_x =   \\
  w_{ ½,0,0}( \x_{1,¼,½} ) \  v^x_{ ½,0,0} +  \\
  w_{1½,0,0}( \x_{1,¼,½} ) \  v^x_{1½,0,0} +  \\
  \vdots \\
  w_{1½,1,1}( \x_{1,¼,½} )\ v^x_{1½,1,1}.
\\]

Since we need to do these for the $x$-component of each input normal, we will
assemble a sparse matrix $\W^x ∈ n × (n_x-1)n_yn_z$ that _interpolates_
$\v^x$ at each point $\p$:

\\[
  ( \W^x \v^x ) ∈ \R^{n×1}
\\]

the transpose of $\W^x$ is not quite its
[_inverse_](https://en.wikipedia.org/wiki/Invertible_matrix), but instead can
be interpreted as _distributing_ values onto staggered grid locations where
$\v^x$ lives:

\\[
  \v^x = (\W^x)^\transpose \N^x.
\\]

Using this definition of $\v^x$ and analogously for $\v^y$ and $\v^z$ we can
construct the vector $\v$ in our energy minimization problem above.

> ### BTW, what's [Poisson](https://en.wikipedia.org/wiki/Siméon_Denis_Poisson) got to do with it?
> 
> The discrete energy minimization problem we've written looks like the squared
> norm of some gradients. An analogous energy in the smooth world is the
> [Dirichlet energy](https://en.wikipedia.org/wiki/Dirichlet's_energy):
>
> \\[
>   E(g) = ∫_Ω ‖∇g‖² dA
> \\]
>
> to _minimize_ this energy with respect to $g$ as an unknown _function_, we
> need to invoke [Calculus of
> Variations](https://en.wikipedia.org/wiki/Calculus_of_variations) and
> [Green's First
> Identity](https://en.wikipedia.org/wiki/Green's_identities#Green.27s_first_identity).
> In doing so we find that minimizers will satisfy:
>
> \\[
>   ∇⋅∇ g = 0 \text{ on Ω},
> \\]
> 
> known as [Laplaces'
> Equation](https://en.wikipedia.org/wiki/Laplace%27s_equation).
>
> If we instead start with a slightly different energy:
> 
> \\[
>   E(g) = ∫_Ω ‖∇g - V‖² dA,
> \\]
> 
> where $V$ is a vector-valued function. Then applying the same machinery we
> find that minimizers will satisfy:
>
> \\[
>   ∇⋅∇ g = ∇⋅V \text{ on Ω},
> \\]
>
> known as [Poisson's
> equation](https://en.wikipedia.org/wiki/Poisson%27s_equation). 
>
>
> Notice that if we interpret the transpose of our gradient matrix
> $\G^\transpose$ as a _divergence matrix_ (we can and we should), then the
> structure of these smooth energies and equations are directly preserved in
> our discrete energies and equations.
>
> This kind of _structure preservation_ is a major criterion for judging
> discrete methods.
>

### Choosing a good iso-value

Constant functions have no gradient. This means that we can add a constant
function to our implicit function $g$ without changing its gradient:

\\[
∇g = ∇(g+c) = ∇g + ∇c = ∇g + 0.
\\]

The same is true for our discrete gradient matrix $\G$: if the vector of grid
values $\g$ is constant then $\G \g$ will be a vector zeros.

This is potentially problematic for our least squares solve: there are many
solutions, since we can just add a constant. Fortunately, we _don't really
care_. It's elegant to say that our surface is defined at $g=0$, but we'd be
just as happy declaring that our surface is defined at $g=c$.

To this end we just need to find _a solution_ $\g$, and then to pick a good
iso-value $σ$.

As suggested in \[Kazhdan et al. 2006\], we can pick a good iso-value by
interpolating our solution $\g$ at each of the input points (since we know
they're on the surface) and averaging their values. In matrix form using our
interpolation matrix $\W$:

\\[
σ = \frac{1}{n} \mathbf{1}^\transpose \W \g,
\\]

where $\mathbf{1} ∈ \R^{n×1}$ is a vector of ones.

### Just how simplified from \[Kazhdan et al. 2006\] is this assignment

Besides the insights above, a major contribution of \[Kazhdan et al. 2006\] was
to setup and solve this problem on an [adaptive
grid](https://en.wikipedia.org/wiki/Adaptive_mesh_refinement) rather than a
regular grid. They also track "confidence" of their input data effecting how
they smooth and interpolate values. As a result, their method is one of the
most highly used surface reconstruction techniques to this day.

> Consider adding your own insights to the wikipedia entry for [this
> method](https://en.wikipedia.org/wiki/Poisson's%5Fequation#Surface%5Freconstruction).

## Tasks

### Read \[Kazhdan et al. 2006\]

This reading task is not directly graded, but it's expected that you read and
understand this paper before moving on to the other tasks.

### `src/fd_interpolate.cpp`

Given a regular finite-difference grid described by the number of nodes on each
side (`nx`, `ny` and `nz`), the grid spacing (`h`), and the location of the
bottom-left-front-most corner node (`corner`), and a list of points (`P`),
construct a sparse matrix `W` of trilinear interpolation weights so that `P = W
* x`. 

### `src/fd_partial_derivative.cpp`

Given a regular finite-difference grid described by the number of nodes on each
side (`nx`, `ny` and `nz`), the grid spacing (`h`), and a desired direction,
construct a sparse matrix `D` to compute first partial derivatives in the given
direction onto the _staggered grid_ in that direction.

### `src/fd_grad.cpp`

Given a regular finite-difference grid described by the number of nodes on each
side (`nx`, `ny` and `nz`), and the grid spacing (`h`), construct a sparse
matrix `G` to compute gradients with each component on its respective staggered
grid.

**_Hint:_** use `fd_partial_derivative` to compute `Dx`, `Dy`, and `Dz` and
then simply concatenate these together to make `G`.

### `src/poisson_surface_reconstruction.cpp`

Given a list of points `P` and the list of corresponding normals `N`, construct
a implicit function `g` over a regular grid (_built for you_) using approach
described above.

You will need to _distribute_ the given normals `N` onto the staggered grid
values in `v` via sparse trilinear interpolation matrices `Wx`, `Wy` and `Wz`
for each staggered grid. 

Then you will need to construct and solve the linear system $\G^\transpose \G
\g = \G^\transpose v$.

Determine the iso-level `sigma` to extract from the `g`.

Feed this implicit function `g` to `igl::copyleft::marching_cubes` to contour
this function into a triangle mesh `V` and `F`.

Make use of `fd_interpolate` and `fd_grad`.

**_Hint:_** Eigen has many different [sparse matrix
solvers](https://eigen.tuxfamily.org/dox-devel/group__TopicSparseSystems.html).
For these _very regular_ matrices, it seems that the [conjugate gradient
method](https://en.wikipedia.org/wiki/Conjugate_gradient_method) will
outperform direct methods such as [Cholesky
factorization](https://en.wikipedia.org/wiki/Cholesky_decomposition). Try
`Eigen::BiCGSTAB`.

**_Hint:_** Debug in debug mode with assertions enabled. For Unix users on the
command line use: 

    cmake -DCMAKE_BUILD_TYPE=Debug ../

but then try out your code in _release mode_ for much better performance

    cmake -DCMAKE_BUILD_TYPE=Release ../
