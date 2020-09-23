# Geometry Processing – Mesh Reconstruction

> **To get started:** Clone this repository then issue
> 
>     git clone --recursive http://github.com/alecjacobson/geometry-processing-mesh-reconstruction.git
>

## Installation, Layout, and Compilation

See
[introduction](http://github.com/alecjacobson/geometry-processing-introduction).

## Execution

Once built, you can execute the assignment from inside the `build/` using 

    ./mesh-reconstruction [path to point cloud]

## Background

In this assignment, we will be implementing a simplified version of the method
in  "Poisson Surface Reconstruction" by Kazhdan et al. 2006. (Your first "task"
will be to read and understand this paper).

Many scanning technologies output a set of <img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/> point samples <img src="/tex/384591906555413c452c93e493b2d4ec.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92230829999999pt height=22.55708729999998pt/> on the 
surface of the object in question. From these points and perhaps the location
of the camera, one can also estimate normals <img src="/tex/bccab73005d96290c8ef588703533a21.svg?invert_in_darkmode&sanitize=true" align=middle width=14.794451099999991pt height=22.55708729999998pt/> to the surface for each point
<img src="/tex/f4e5add43718dd0d01e8a2cf425322f8.svg?invert_in_darkmode&sanitize=true" align=middle width=43.515672749999986pt height=22.55708729999998pt/>. This image shows the `data/elephant.pwn` input
data with a white dot for each point and a red line segment pointing outward for
each corresponding normal vector.

![](images/elephant-points-normals.jpg)

For shape analysis, visualization and other downstream geometry processing
phases, we would like to convert this finitely sampled _point cloud_ data into
an _explicit continuous surface representation_: i.e., a [triangle
mesh](https://en.wikipedia.org/wiki/Triangulation%5F(topology)) (a special case
of a [polygon mesh](https://en.wikipedia.org/wiki/Polygon%5Fmesh)). This image
shows the corresponding output mesh for `data/elephant.pwn` input data above:

![](images/elephant-mesh.jpg)

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
[level-set](https://en.wikipedia.org/wiki/Level%5Fset) of some function <img src="/tex/bb9905f60ed7973f6bb9b104a3cfb901.svg?invert_in_darkmode&sanitize=true" align=middle width=83.42985254999998pt height=26.76175259999998pt/> mapping all points in space to a scalar value. For example, we may define
the surface <img src="/tex/cffca7c1479a119f8929bf5726372e56.svg?invert_in_darkmode&sanitize=true" align=middle width=20.142632399999986pt height=22.831056599999986pt/> of some solid, volumetric shape <img src="/tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode&sanitize=true" align=middle width=10.502226899999991pt height=22.55708729999998pt/> to be all points <img src="/tex/894fe1065bc022997ddadbf01f0cf213.svg?invert_in_darkmode&sanitize=true" align=middle width=50.798807399999994pt height=26.76175259999998pt/> such that <img src="/tex/7692e3f86410f1343c20d3bc0f95b02c.svg?invert_in_darkmode&sanitize=true" align=middle width=62.51131094999999pt height=24.65753399999998pt/>, where we may arbitrarily set <img src="/tex/5439c1794588f864d0fcccab4d28f058.svg?invert_in_darkmode&sanitize=true" align=middle width=40.42565669999999pt height=27.77565449999998pt/>.

<p align="center"><img src="/tex/437fbd49ccec7ed33a81f87d3df7d43a.svg?invert_in_darkmode&sanitize=true" align=middle width=182.34527025pt height=18.312383099999998pt/></p>


On the computer, it is straightforward
[discretize](https://en.wikipedia.org/wiki/Discretization) an implicit
function. We define a regular 3D grid of
[voxels](http://en.wikipedia.org/wiki/Voxel) containing at least the [bounding
box](https://en.wikipedia.org/wiki/Minimum%5Fbounding%5Fbox) of <img src="/tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode&sanitize=true" align=middle width=10.502226899999991pt height=22.55708729999998pt/>. At each
node in the grid <img src="/tex/4d62a3e4e599acad7cc1af2f73ad0543.svg?invert_in_darkmode&sanitize=true" align=middle width=35.09903264999999pt height=14.611878600000017pt/> we store the value of the implicit function
<img src="/tex/e611a18d5ae8daa7807f07b9889155d0.svg?invert_in_darkmode&sanitize=true" align=middle width=57.13667024999998pt height=24.65753399999998pt/>. This defines <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> _everywhere_ in the grid via [trilinear
interpolation](https://en.wikipedia.org/wiki/Trilinear_interpolation). 

![](images/trilinear-interpolation.jpg)

For example, consider a point <img src="/tex/e89a81a146afd2d4c13849f2ce122376.svg?invert_in_darkmode&sanitize=true" align=middle width=90.7851582pt height=27.77565449999998pt/> lying in the middle of the
bottom-most, front-most, left-most cell. We know the values at the eight
corners. Trilinear interpolation can be understood as [linear
interpolation](https://en.wikipedia.org/wiki/Linear_interpolation) in the
<img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-direction by <img src="/tex/2bbf4eee68fd087c4e89458ffc858ecb.svg?invert_in_darkmode&sanitize=true" align=middle width=6.552545999999997pt height=27.77565449999998pt/> on each <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-axis-aligned edge, resulting in four values
_living_ on the same plane. These can then be linearly interpolated in the <img src="/tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode&sanitize=true" align=middle width=8.649225749999989pt height=14.15524440000002pt/>
direction by <img src="/tex/b1467b64271e0e80be78cc8537d840cb.svg?invert_in_darkmode&sanitize=true" align=middle width=6.552545999999997pt height=27.77565449999998pt/> resulting in two points on the same line, and finally in the
<img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/> direction by <img src="/tex/2eacd907c22e9ad1db070a00eef0f2cd.svg?invert_in_darkmode&sanitize=true" align=middle width=6.552545999999997pt height=27.77565449999998pt/> to get to our evaluation point <img src="/tex/34443414f341f1768e4a3663859ea298.svg?invert_in_darkmode&sanitize=true" align=middle width=58.890412349999984pt height=27.77565449999998pt/>.

An implicit surface stored as the level-set of a trilinearly interpolated grid
can be _contoured_ into a triangle mesh via the [Marching Cubes
Algorithm](https://en.wikipedia.org/wiki/Marching&5Fcubes). 
For the purposes of this assignment, we will treat this as a [black
box](https://en.wikipedia.org/wiki/Black%5Fbox). Instead, we focus on
determining what values for <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> to store on the grid.

## Characteristic functions of solids

We assume that our set of points <img src="/tex/384591906555413c452c93e493b2d4ec.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92230829999999pt height=22.55708729999998pt/> lie on the surface <img src="/tex/cffca7c1479a119f8929bf5726372e56.svg?invert_in_darkmode&sanitize=true" align=middle width=20.142632399999986pt height=22.831056599999986pt/> of some physical
[solid](https://en.wikipedia.org/wiki/Solid) object <img src="/tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode&sanitize=true" align=middle width=10.502226899999991pt height=22.55708729999998pt/>. This solid object
must have some non-trivial volume that we can calculate _abstractly_ as the
integral of unit density over the solid:

<p align="center"><img src="/tex/17fd0995c494e4eb8497758396409438.svg?invert_in_darkmode&sanitize=true" align=middle width=57.414333899999995pt height=47.164758299999995pt/></p>


We can rewrite this definite integral as an indefinite integral over all of
<img src="/tex/1e500ed0738135b9e26f388580e20b32.svg?invert_in_darkmode&sanitize=true" align=middle width=20.730553799999992pt height=26.76175259999998pt/>:

<p align="center"><img src="/tex/d704526aa8d189d72fa711747b30ed72.svg?invert_in_darkmode&sanitize=true" align=middle width=96.0056361pt height=47.85421245pt/></p>


by introducing the [characteristic
function](https://en.wikipedia.org/wiki/Indicator%5Ffunction) of <img src="/tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode&sanitize=true" align=middle width=10.502226899999991pt height=22.55708729999998pt/>, that is
_one_ for points inside of the shape and _zero_ for points outside of <img src="/tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode&sanitize=true" align=middle width=10.502226899999991pt height=22.55708729999998pt/>:

<p align="center"><img src="/tex/fc0332b4de4ffb25b3d11d449336566b.svg?invert_in_darkmode&sanitize=true" align=middle width=179.71458779999998pt height=49.315569599999996pt/></p>


Compared to typical [implicit surface
functions](https://en.wikipedia.org/wiki/Implicit%5Fsurface), this function
represents the surface <img src="/tex/cffca7c1479a119f8929bf5726372e56.svg?invert_in_darkmode&sanitize=true" align=middle width=20.142632399999986pt height=22.831056599999986pt/> of the shape <img src="/tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode&sanitize=true" align=middle width=10.502226899999991pt height=22.55708729999998pt/> as the _discontinuity_ between
the one values and the zero values. Awkwardly, the gradient of the
characteristic function <img src="/tex/8fc2223d19b0886790872904d6b97578.svg?invert_in_darkmode&sanitize=true" align=middle width=32.24885729999999pt height=22.465723500000017pt/> is _not defined_ along <img src="/tex/cffca7c1479a119f8929bf5726372e56.svg?invert_in_darkmode&sanitize=true" align=middle width=20.142632399999986pt height=22.831056599999986pt/>.

One of the key observations made in [Kazhdan et al. 2006] is that the gradient
of a infinitesimally [mollified](https://en.wikipedia.org/wiki/Mollifier)
(smoothed) characteristic function: 

  1. points in the direction of the normal near the surface <img src="/tex/cffca7c1479a119f8929bf5726372e56.svg?invert_in_darkmode&sanitize=true" align=middle width=20.142632399999986pt height=22.831056599999986pt/>, and 
  2. is zero everywhere else.

Our goal will be to use our points <img src="/tex/384591906555413c452c93e493b2d4ec.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92230829999999pt height=22.55708729999998pt/> and normals <img src="/tex/bccab73005d96290c8ef588703533a21.svg?invert_in_darkmode&sanitize=true" align=middle width=14.794451099999991pt height=22.55708729999998pt/> to _optimize_ an
implicit function <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> over a regular grid, so that its gradient <img src="/tex/ecb02abb9d3c72a23e6c3e776b8f3b53.svg?invert_in_darkmode&sanitize=true" align=middle width=22.129049249999987pt height=22.465723500000017pt/> meets
these two criteria. In that way, our <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> will be an approximation of the
mollified characteristic function.

## Poisson surface reconstruction

### Or: how I learned to stop worrying and minimize squared Gradients

<span id=assumptions>
Let us start by making two assumptions:
</span>

  1. we know how to compute <img src="/tex/ecb02abb9d3c72a23e6c3e776b8f3b53.svg?invert_in_darkmode&sanitize=true" align=middle width=22.129049249999987pt height=22.465723500000017pt/> at each node location <img src="/tex/4d62a3e4e599acad7cc1af2f73ad0543.svg?invert_in_darkmode&sanitize=true" align=middle width=35.09903264999999pt height=14.611878600000017pt/>, and
  2. our input points <img src="/tex/384591906555413c452c93e493b2d4ec.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92230829999999pt height=22.55708729999998pt/> all lie perfectly at grid nodes: 
  <img src="/tex/ffd3200af88dde1e39d8a5b6c2cc16ec.svg?invert_in_darkmode&sanitize=true" align=middle width=88.45497209999998pt height=22.831056599999986pt/>.

We will find out these assumptions are not realistic and we will have to relax
them (i.e., we **_will not_** make these assumptions in the completion of the
tasks). However, it will make the following algorithmic description easier on
the first pass.

If our points <img src="/tex/384591906555413c452c93e493b2d4ec.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92230829999999pt height=22.55708729999998pt/> lie at grid points, then our corresponding normals <img src="/tex/bccab73005d96290c8ef588703533a21.svg?invert_in_darkmode&sanitize=true" align=middle width=14.794451099999991pt height=22.55708729999998pt/> also
_live_ at grid points. This leads to a very simple set of linear equations to
define a function <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> with a gradient equal to the surface normal at the
surface and zero gradient away from the surface:

<p align="center"><img src="/tex/8c063c3e6477d6dc92852a303442db8f.svg?invert_in_darkmode&sanitize=true" align=middle width=354.0237888pt height=139.06941179999998pt/></p>


This is a _vector-valued_ equation. The gradients, normals and zero-vectors are
three-dimensional (e.g., <img src="/tex/97264e77606b79af66ff5c3064a3d179.svg?invert_in_darkmode&sanitize=true" align=middle width=62.950720799999985pt height=26.76175259999998pt/>). In effect, this is _three equations_ for
every grid node.

Since we only need a single number at each grid node (the value of <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/>), we
have _too many_ equations.

Like many geometry processing algorithms confronted with such an [over
determined](https://en.wikipedia.org/wiki/Overdetermined%5Fsystem), we will
_optimize_ for the solution that best _minimizes_ the error of equation:

<p align="center"><img src="/tex/8307d3640b2448d0bee6f177cbda22dc.svg?invert_in_darkmode&sanitize=true" align=middle width=155.22651374999998pt height=18.905967299999997pt/></p>


We will treat the error of each grid location equally by minimizing the sum
over all grid locations:

<p align="center"><img src="/tex/f5c9eb9c5fad33ab1f8895e3863d1d88.svg?invert_in_darkmode&sanitize=true" align=middle width=276.97991145pt height=43.346758949999995pt/></p>


where <img src="/tex/5f3cc59831e6b4aef298a2dacada3fe7.svg?invert_in_darkmode&sanitize=true" align=middle width=9.714576299999992pt height=14.611878600000017pt/> (written in boldface) is a vector of _unknown_ grid-nodes values,
where <img src="/tex/c4ca9607e0256b282b5cdabc3a37e938.svg?invert_in_darkmode&sanitize=true" align=middle width=112.83864404999998pt height=24.65753399999998pt/>. 

Part of the convenience of working on a regular grid is that we can use the
[finite difference
method](https://en.wikipedia.org/wiki/Finite_difference_method) to approximate
the gradient <img src="/tex/ecb02abb9d3c72a23e6c3e776b8f3b53.svg?invert_in_darkmode&sanitize=true" align=middle width=22.129049249999987pt height=22.465723500000017pt/> on the grid.

After revisiting [our assumptions](#assumptions), we will be able to compute
approximations of 
the <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-, <img src="/tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode&sanitize=true" align=middle width=8.649225749999989pt height=14.15524440000002pt/>- and <img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/>-components of <img src="/tex/ecb02abb9d3c72a23e6c3e776b8f3b53.svg?invert_in_darkmode&sanitize=true" align=middle width=22.129049249999987pt height=22.465723500000017pt/> via a [sparse
matrix](https://en.wikipedia.org/wiki/Sparse%5Fmatrix) multiplication of
a "gradient matrix" <img src="/tex/ae44fa1818647ed39ba79b316ce5dd87.svg?invert_in_darkmode&sanitize=true" align=middle width=14.86294424999999pt height=22.55708729999998pt/> and our vector of unknown grid values <img src="/tex/5f3cc59831e6b4aef298a2dacada3fe7.svg?invert_in_darkmode&sanitize=true" align=middle width=9.714576299999992pt height=14.611878600000017pt/>. We will be
able to write the
minimization problem above in matrix form:

<p align="center"><img src="/tex/eea7dd94c059c1b0a4871b27e5226110.svg?invert_in_darkmode&sanitize=true" align=middle width=125.58889695pt height=33.814738649999995pt/></p>


or equivalently after expanding the norm:

<p align="center"><img src="/tex/764485b6489c166576a9c344e1f9cec2.svg?invert_in_darkmode&sanitize=true" align=middle width=268.53430725pt height=33.814738649999995pt/></p>


This is a quadratic "energy" function of the variables of <img src="/tex/5f3cc59831e6b4aef298a2dacada3fe7.svg?invert_in_darkmode&sanitize=true" align=middle width=9.714576299999992pt height=14.611878600000017pt/>, its minimum occurs when
an infinitesimal change in <img src="/tex/5f3cc59831e6b4aef298a2dacada3fe7.svg?invert_in_darkmode&sanitize=true" align=middle width=9.714576299999992pt height=14.611878600000017pt/> produces no change in the energy:

<p align="center"><img src="/tex/3ae9f01dc94f6f227393748350a53b77.svg?invert_in_darkmode&sanitize=true" align=middle width=208.9484826pt height=37.0084374pt/></p>


Applying this derivative gives us a _sparse_ system of linear equations

<p align="center"><img src="/tex/de359c205c9bfe1240b157e694922b0f.svg?invert_in_darkmode&sanitize=true" align=middle width=109.3719528pt height=17.9744895pt/></p>


We will assume that we can solve this using a black box sparse solver.

Now, let's revisit [our assumptions](#assumptions).

### Gradients on a regular grid

The gradient of a function <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> in 3D is nothing more than a vector containing
partial derivatives in each coordinate direction:

<p align="center"><img src="/tex/492a25e9363c8fec347ccd7b6ed47b5f.svg?invert_in_darkmode&sanitize=true" align=middle width=155.9496675pt height=69.3006039pt/></p>


We will approximate each partial derivative individually. Let's consider the
partial derivative in the <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/> direction, <img src="/tex/eab7528ff6cda4ba4cabe93ba6c4d2ea.svg?invert_in_darkmode&sanitize=true" align=middle width=68.08791659999999pt height=24.65753399999998pt/>, and we will assume
without loss of generality that what we derive applies _symmetrically_ for <img src="/tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode&sanitize=true" align=middle width=8.649225749999989pt height=14.15524440000002pt/>
and <img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/>.

The partial derivative in the <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-direction is a one-dimensional derivative.
This couldn't be easier to do with finite differences. We approximate the
derivative of the function <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> with respect to the <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/> direction is the
difference between the function evaluated at one grid node and at the grid node
_before_ it in the <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-direction then divided by the spatial distance between
adjacent nodes <img src="/tex/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode&sanitize=true" align=middle width=9.47111549999999pt height=22.831056599999986pt/> (i.e., the grid step size):

<p align="center"><img src="/tex/47486f504760931f80d59be4c7e5f4d4.svg?invert_in_darkmode&sanitize=true" align=middle width=223.47799815pt height=38.02854825pt/></p>


where we use the index <img src="/tex/bd8a3a3788c7e885e5bb709066e6888a.svg?invert_in_darkmode&sanitize=true" align=middle width=34.27956179999999pt height=27.77565449999998pt/> to indicate that this derivative in the
<img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-direction lives on a [staggered
grid](https://en.wikipedia.org/wiki/Staggered%5Fgrid) _in between_ the grid
nodes where the function values for <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/>.

The following pictures show a 2D example, where <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> lives on the nodes of a
<img src="/tex/3e2a76e23fe9c67f0c6a1d099d494151.svg?invert_in_darkmode&sanitize=true" align=middle width=36.52961069999999pt height=21.18721440000001pt/> blue grid:

![](images/primary-grid.jpg)

The partial derivatives of <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> with respect to the <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-direction <img src="/tex/eab7528ff6cda4ba4cabe93ba6c4d2ea.svg?invert_in_darkmode&sanitize=true" align=middle width=68.08791659999999pt height=24.65753399999998pt/>
live on a <img src="/tex/6e9103d18236e9931449ce1f3a303f5c.svg?invert_in_darkmode&sanitize=true" align=middle width=36.52961069999999pt height=21.18721440000001pt/> green, staggered grid:

![](images/staggered-grid-x.jpg)

The partial derivatives of <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> with respect to the <img src="/tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode&sanitize=true" align=middle width=8.649225749999989pt height=14.15524440000002pt/>-direction <img src="/tex/3236f157a53b70469c7a6de44c7d21f0.svg?invert_in_darkmode&sanitize=true" align=middle width=67.34215289999999pt height=24.65753399999998pt/>
live on a <img src="/tex/0ee0c8e0c23693d9c6b800197f53ea56.svg?invert_in_darkmode&sanitize=true" align=middle width=36.52961069999999pt height=21.18721440000001pt/> yellow, staggered grid:

![](images/staggered-grid-x-and-y.jpg)

Letting <img src="/tex/793a1b090156858c6786a258487ce31c.svg?invert_in_darkmode&sanitize=true" align=middle width=106.42512869999997pt height=26.76175259999998pt/> be column vector of function values on the
_primary grid_ (blue in the example pictures), we can construct a sparse matrix
<img src="/tex/54598fd8e0d7e440dfa0c17e17f8be10.svg?invert_in_darkmode&sanitize=true" align=middle width=184.82548755pt height=29.190975000000005pt/> so that each row <img src="/tex/4772bc64bc27c6be242900b098354180.svg?invert_in_darkmode&sanitize=true" align=middle width=156.14297655pt height=26.76175259999998pt/> computes the partial derivative at the corresponding
staggered grid location <img src="/tex/73adb055289ee2e3d0fed569b21285ff.svg?invert_in_darkmode&sanitize=true" align=middle width=54.911945549999984pt height=16.026044099999982pt/>. The <img src="/tex/d30a65b936d8007addc9c789d5a7ae49.svg?invert_in_darkmode&sanitize=true" align=middle width=6.849367799999992pt height=22.831056599999986pt/>th entry in that row receives a
value only for neighboring primary grid nodes:

<p align="center"><img src="/tex/dcf8770cd411ce907b99ca3bea200f9e.svg?invert_in_darkmode&sanitize=true" align=middle width=251.90811029999998pt height=69.0417981pt/></p>


> #### Indexing 3D arrays
> 
> Now, obviously in our code we cannot _index_ the column vector <img src="/tex/5f3cc59831e6b4aef298a2dacada3fe7.svg?invert_in_darkmode&sanitize=true" align=middle width=9.714576299999992pt height=14.611878600000017pt/> by a
> triplet of numbers <img src="/tex/daa72bcc65930e39126b00858df4e55c.svg?invert_in_darkmode&sanitize=true" align=middle width=52.58594714999999pt height=24.65753399999998pt/> or the rows of <img src="/tex/98a6924bf931c3807f8deddadd870393.svg?invert_in_darkmode&sanitize=true" align=middle width=21.95201414999999pt height=22.55708729999998pt/> by the triplet
> <img src="/tex/b999fa3230bc59907de3384d2a0a1ab2.svg?invert_in_darkmode&sanitize=true" align=middle width=66.73646264999999pt height=27.77565449999998pt/>. We will assume that <img src="/tex/8a567299f0c8eee3484004a8814e0a9d.svg?invert_in_darkmode&sanitize=true" align=middle width=34.573921799999994pt height=14.611878600000017pt/> refers to
> `g(i+j*n_x+k*n_y*n_x)`. Similarly, for the staggered grid subscripts
> <img src="/tex/b999fa3230bc59907de3384d2a0a1ab2.svg?invert_in_darkmode&sanitize=true" align=middle width=66.73646264999999pt height=27.77565449999998pt/> we will assume that <img src="/tex/d2b4c9acbd841c2d920f43cbd427168b.svg?invert_in_darkmode&sanitize=true" align=middle width=79.88909444999999pt height=24.65753399999998pt/> refers to the matrix
> entry `Dx(i+j*n_x+k*n_y*n_x,l)`, where the <img src="/tex/bd8a3a3788c7e885e5bb709066e6888a.svg?invert_in_darkmode&sanitize=true" align=middle width=34.27956179999999pt height=27.77565449999998pt/> has been _rounded down_.
>

We can similarly build matrices <img src="/tex/81ee55d5839b20e95d03b438b3a60dbd.svg?invert_in_darkmode&sanitize=true" align=middle width=21.57724964999999pt height=22.55708729999998pt/> and <img src="/tex/2f11c6120fcb1f7a8dd90b5a25a7216b.svg?invert_in_darkmode&sanitize=true" align=middle width=21.250000199999988pt height=22.55708729999998pt/> and _stack_ these matrices
vertically to create a gradient matrix <img src="/tex/ae44fa1818647ed39ba79b316ce5dd87.svg?invert_in_darkmode&sanitize=true" align=middle width=14.86294424999999pt height=22.55708729999998pt/>:

<p align="center"><img src="/tex/3cec78ebec6c7208cc0c77bdbe913bf8.svg?invert_in_darkmode&sanitize=true" align=middle width=442.699092pt height=59.1786591pt/></p>


This implies that our vector <img src="/tex/f6fc3ac36dff143d4aac9d145fadc77e.svg?invert_in_darkmode&sanitize=true" align=middle width=10.239687149999991pt height=14.611878600000017pt/> of zeros and normals in our minimization
problem should not _live_ on the primary, but rather it, too, should be broken
into <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-, <img src="/tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode&sanitize=true" align=middle width=8.649225749999989pt height=14.15524440000002pt/>- and <img src="/tex/f93ce33e511096ed626b4719d50f17d2.svg?invert_in_darkmode&sanitize=true" align=middle width=8.367621899999993pt height=14.15524440000002pt/>-components that live of their resepctive staggered
grids:

<p align="center"><img src="/tex/517cbc8a282a06299b941df2bf4ad9f2.svg?invert_in_darkmode&sanitize=true" align=middle width=400.9655694pt height=59.1786591pt/></p>


This leads to addressing our second assumption.

### B-b-b-b-but the input normals might not be at grid node locations?

At this point, we would _actually_ liked to have had that our input normals
were given component-wise on the staggered grid. Then we could immediate stick
them into <img src="/tex/f6fc3ac36dff143d4aac9d145fadc77e.svg?invert_in_darkmode&sanitize=true" align=middle width=10.239687149999991pt height=14.611878600000017pt/>. But this doesn't make much sense as each normal <img src="/tex/985e9d424986b945b21cbfc30f69ee2c.svg?invert_in_darkmode&sanitize=true" align=middle width=16.00457264999999pt height=14.611878600000017pt/> _lives_
at its associated point <img src="/tex/abba53cfa4754e0fb3447fe5e0be0b78.svg?invert_in_darkmode&sanitize=true" align=middle width=16.00457264999999pt height=14.611878600000017pt/>, regardless of any grids.

To remedy this, we will distribute each component of each input normal <img src="/tex/985e9d424986b945b21cbfc30f69ee2c.svg?invert_in_darkmode&sanitize=true" align=middle width=16.00457264999999pt height=14.611878600000017pt/>
to <img src="/tex/f6fc3ac36dff143d4aac9d145fadc77e.svg?invert_in_darkmode&sanitize=true" align=middle width=10.239687149999991pt height=14.611878600000017pt/> at the corresponding staggered grid node location.

For example, consider the normal <img src="/tex/b56595d2a30a0af329086562ca12d521.svg?invert_in_darkmode&sanitize=true" align=middle width=10.502226899999991pt height=14.611878600000017pt/> at some point <img src="/tex/f8418a00a5940164cc6a57d88beceeab.svg?invert_in_darkmode&sanitize=true" align=middle width=41.443117649999984pt height=16.026044099999982pt/>. Conceptually,
we'll think of the <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-component of the normal <img src="/tex/322d8f61a96f4dd07a0c599482268dfe.svg?invert_in_darkmode&sanitize=true" align=middle width=17.32124954999999pt height=14.15524440000002pt/> as floating in the
staggered grid corresponding to <img src="/tex/98a6924bf931c3807f8deddadd870393.svg?invert_in_darkmode&sanitize=true" align=middle width=21.95201414999999pt height=22.55708729999998pt/>, in between the eight staggered grid
locations:

<p align="center"><img src="/tex/11420dc7bbad910d2c10985d1b742ea8.svg?invert_in_darkmode&sanitize=true" align=middle width=443.98231185000003pt height=19.42855035pt/></p>


Each of these staggered grid nodes has a corresponding <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/> value in the vector
<img src="/tex/2721b133f6f6619c279972b564b09689.svg?invert_in_darkmode&sanitize=true" align=middle width=17.69405714999999pt height=21.839370299999988pt/>.

We will distribute <img src="/tex/322d8f61a96f4dd07a0c599482268dfe.svg?invert_in_darkmode&sanitize=true" align=middle width=17.32124954999999pt height=14.15524440000002pt/> to these entries in <img src="/tex/2721b133f6f6619c279972b564b09689.svg?invert_in_darkmode&sanitize=true" align=middle width=17.69405714999999pt height=21.839370299999988pt/> by _adding_ a partial amount
of <img src="/tex/322d8f61a96f4dd07a0c599482268dfe.svg?invert_in_darkmode&sanitize=true" align=middle width=17.32124954999999pt height=14.15524440000002pt/> to each. I.e., 

<p align="center"><img src="/tex/2ca7f2cd2542e5d50a1ca0023d729ff8.svg?invert_in_darkmode&sanitize=true" align=middle width=641.4843055499999pt height=29.58934275pt/></p>

where <img src="/tex/76c0bc21c84b795247769888ea687c90.svg?invert_in_darkmode&sanitize=true" align=middle width=80.63019194999998pt height=24.65753399999998pt/> is the trilinear interpolation _weight_ associate with
staggered grid node <img src="/tex/f5635624df56862986d167f2d91cc6e6.svg?invert_in_darkmode&sanitize=true" align=middle width=54.72929879999998pt height=16.026044099999982pt/> to interpolate a value at the point <img src="/tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg?invert_in_darkmode&sanitize=true" align=middle width=10.502226899999991pt height=14.611878600000017pt/>.
The trilinear interpolation weights so that:

<p align="center"><img src="/tex/8ff5d71a18238a4be365916ca7d05d68.svg?invert_in_darkmode&sanitize=true" align=middle width=549.94754265pt height=24.403666649999998pt/></p>


Since we need to do these for the <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-component of each input normal, we will
assemble a sparse matrix <img src="/tex/da5edae1fb4f79cabe92a9ab16962ad7.svg?invert_in_darkmode&sanitize=true" align=middle width=171.75790995pt height=24.65753399999998pt/> that _interpolates_
<img src="/tex/2721b133f6f6619c279972b564b09689.svg?invert_in_darkmode&sanitize=true" align=middle width=17.69405714999999pt height=21.839370299999988pt/> at each point <img src="/tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg?invert_in_darkmode&sanitize=true" align=middle width=10.502226899999991pt height=14.611878600000017pt/>:

<p align="center"><img src="/tex/7da1a14aa73a72098c04b81d8feed597.svg?invert_in_darkmode&sanitize=true" align=middle width=118.60523565pt height=18.312383099999998pt/></p>


the transpose of <img src="/tex/473d4a136e73dec1b14f253176da2f53.svg?invert_in_darkmode&sanitize=true" align=middle width=27.26022089999999pt height=22.55708729999998pt/> is not quite its
[_inverse_](https://en.wikipedia.org/wiki/Invertible_matrix), but instead can
be interpreted as _distributing_ values onto staggered grid locations where
<img src="/tex/2721b133f6f6619c279972b564b09689.svg?invert_in_darkmode&sanitize=true" align=middle width=17.69405714999999pt height=21.839370299999988pt/> lives:

<p align="center"><img src="/tex/7ff3dcedb9ca37b4fa9ef4b6aeb87111.svg?invert_in_darkmode&sanitize=true" align=middle width=118.11057719999998pt height=18.88772655pt/></p>


Using this definition of <img src="/tex/2721b133f6f6619c279972b564b09689.svg?invert_in_darkmode&sanitize=true" align=middle width=17.69405714999999pt height=21.839370299999988pt/> and analogously for <img src="/tex/05af85763c308fe9dfd93f44e8149a59.svg?invert_in_darkmode&sanitize=true" align=middle width=17.31929099999999pt height=21.839370299999988pt/> and <img src="/tex/41bf049c525434295c5c91ab03d41324.svg?invert_in_darkmode&sanitize=true" align=middle width=16.99204154999999pt height=21.839370299999988pt/> we can
construct the vector <img src="/tex/f6fc3ac36dff143d4aac9d145fadc77e.svg?invert_in_darkmode&sanitize=true" align=middle width=10.239687149999991pt height=14.611878600000017pt/> in our energy minimization problem above.

> ### BTW, what's [Poisson](https://en.wikipedia.org/wiki/Siméon_Denis_Poisson) got to do with it?
> 
> The discrete energy minimization problem we've written looks like the squared
> norm of some gradients. An analogous energy in the smooth world is the
> [Dirichlet energy](https://en.wikipedia.org/wiki/Dirichlet's_energy):
>
> <p align="center"><img src="/tex/99d0a596e92be38901109a081d07d64f.svg?invert_in_darkmode&sanitize=true" align=middle width=145.11989745pt height=37.3519608pt/></p>
>
>
> to _minimize_ this energy with respect to <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> as an unknown _function_, we
> need to invoke [Calculus of
> Variations](https://en.wikipedia.org/wiki/Calculus_of_variations) and
> [Green's First
> Identity](https://en.wikipedia.org/wiki/Green's_identities#Green.27s_first_identity).
> In doing so we find that minimizers will satisfy:
>
> <p align="center"><img src="/tex/0c8c492bada7dbc998c2939499c8dd3f.svg?invert_in_darkmode&sanitize=true" align=middle width=122.58544319999999pt height=14.42921205pt/></p>
>
> 
> known as [Laplaces'
> Equation](https://en.wikipedia.org/wiki/Laplace%27s_equation).
>
> If we instead start with a slightly different energy:
> 
> <p align="center"><img src="/tex/66db0c1cae14fe5e2a13f61b3b03d67a.svg?invert_in_darkmode&sanitize=true" align=middle width=183.0193431pt height=37.3519608pt/></p>
>
> 
> where <img src="/tex/a9a3a4a202d80326bda413b5562d5cd1.svg?invert_in_darkmode&sanitize=true" align=middle width=13.242037049999992pt height=22.465723500000017pt/> is a vector-valued function. Then applying the same machinery we
> find that minimizers will satisfy:
>
> <p align="center"><img src="/tex/8c2d642b999ade726977a143fdb6419e.svg?invert_in_darkmode&sanitize=true" align=middle width=153.1789182pt height=14.42921205pt/></p>
>
> known as [Poisson's
> equation](https://en.wikipedia.org/wiki/Poisson%27s_equation). 
>
>
> Notice that if we interpret the transpose of our gradient matrix
> <img src="/tex/f0c74577699a7ca84f8d507f9fd49d62.svg?invert_in_darkmode&sanitize=true" align=middle width=23.21353484999999pt height=27.91243950000002pt/> as a _divergence matrix_ (we can and we should), then the
> structure of these smooth energies and equations are directly preserved in
> our discrete energies and equations.
>
> This kind of _structure preservation_ is a major criterion for judging
> discrete methods.
>

### Choosing a good iso-value

Constant functions have no gradient. This means that we can add a constant
function to our implicit function <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/> without changing its gradient:

<p align="center"><img src="/tex/6a2b92e3660606561092ff62529e4533.svg?invert_in_darkmode&sanitize=true" align=middle width=268.03973129999997pt height=16.438356pt/></p>


The same is true for our discrete gradient matrix <img src="/tex/ae44fa1818647ed39ba79b316ce5dd87.svg?invert_in_darkmode&sanitize=true" align=middle width=14.86294424999999pt height=22.55708729999998pt/>: if the vector of grid
values <img src="/tex/5f3cc59831e6b4aef298a2dacada3fe7.svg?invert_in_darkmode&sanitize=true" align=middle width=9.714576299999992pt height=14.611878600000017pt/> is constant then <img src="/tex/ba3f16549018145cc787dd0e73df1c2b.svg?invert_in_darkmode&sanitize=true" align=middle width=24.57752054999999pt height=22.55708729999998pt/> will be a vector zeros.

This is potentially problematic for our least squares solve: there are many
solutions, since we can just add a constant. Fortunately, we _don't really
care_. It's elegant to say that our surface is defined at <img src="/tex/40ecdbad5a9ebcc6dae0bf86ffeddff6.svg?invert_in_darkmode&sanitize=true" align=middle width=38.567197349999994pt height=21.18721440000001pt/>, but we'd be
just as happy declaring that our surface is defined at <img src="/tex/4695cb3c40f1831d62ab5a07ecb94255.svg?invert_in_darkmode&sanitize=true" align=middle width=37.46179304999998pt height=14.15524440000002pt/>.

To this end we just need to find _a solution_ <img src="/tex/5f3cc59831e6b4aef298a2dacada3fe7.svg?invert_in_darkmode&sanitize=true" align=middle width=9.714576299999992pt height=14.611878600000017pt/>, and then to pick a good
iso-value <img src="/tex/c888b5dd080ec7f68a997efb3ed6b96b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.98290094999999pt height=14.15524440000002pt/>.

As suggested in \[Kazhdan et al. 2006\], we can pick a good iso-value by
interpolating our solution <img src="/tex/5f3cc59831e6b4aef298a2dacada3fe7.svg?invert_in_darkmode&sanitize=true" align=middle width=9.714576299999992pt height=14.611878600000017pt/> at each of the input points (since we know
they're on the surface) and averaging their values. For an appropriate
interpolation matrix <img src="/tex/380c103b60c66d6420ec8923cdc6e6e8.svg?invert_in_darkmode&sanitize=true" align=middle width=19.80585089999999pt height=22.55708729999998pt/> on the _primary (non-staggered) grid_ this can be
written as:

<p align="center"><img src="/tex/c3f6ab15dfc4f9089c52bc6258c8a0b6.svg?invert_in_darkmode&sanitize=true" align=middle width=98.16117134999999pt height=32.990165999999995pt/></p>


where <img src="/tex/03a56bc7e563c16280f1195482f10205.svg?invert_in_darkmode&sanitize=true" align=middle width=68.67374084999999pt height=26.76175259999998pt/> is a vector of ones.

### Just how much does this assignment simplify \[Kazhdan et al. 2006\]?

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
construct a sparse matrix `W` of trilinear interpolation weights so that 
`P = W * x`. 

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

> #### Hint 
> Use `fd_partial_derivative` to compute `Dx`, `Dy`, and `Dz` and
> then simply concatenate these together to make `G`.

### `src/poisson_surface_reconstruction.cpp`

Given a list of points `P` and the list of corresponding normals `N`, construct
a implicit function `g` over a regular grid (_built for you_) using approach
described above.

You will need to _distribute_ the given normals `N` onto the staggered grid
values in `v` via sparse trilinear interpolation matrices `Wx`, `Wy` and `Wz`
for each staggered grid. 

Then you will need to construct and solve the linear system <img src="/tex/190806857015e55cba624fe6a29d20ec.svg?invert_in_darkmode&sanitize=true" align=middle width=104.8057329pt height=27.91243950000002pt/>.

Determine the iso-level `sigma` to extract from the `g`.

Feed this implicit function `g` to `igl::copyleft::marching_cubes` to contour
this function into a triangle mesh `V` and `F`.

Make use of `fd_interpolate` and `fd_grad`.

> #### Hint
> Eigen has many different [sparse matrix
> solvers](https://eigen.tuxfamily.org/dox-devel/group__TopicSparseSystems.html).
> For these _very regular_ matrices, it seems that the [conjugate gradient
> method](https://en.wikipedia.org/wiki/Conjugate_gradient_method) will
> outperform direct methods such as [Cholesky
> factorization](https://en.wikipedia.org/wiki/Cholesky_decomposition). Try
> `Eigen::BiCGSTAB`.

-------------------------------------------------------------------------------

> #### Hint 
> Debug in debug mode with assertions enabled. For Unix users on the
> command line use: 
> 
>     cmake -DCMAKE_BUILD_TYPE=Debug ../
> 
> but then try out your code in _release mode_ for much better performance
> 
>     cmake -DCMAKE_BUILD_TYPE=Release ../
