# Multi-Scale Curvature Estimation

This project implements a novel algorithm for computing multi-scale mean curvature fields of triangle meshes.

For the mean curvature computation at a given vertex and radius, integral invariants with the ball neighborhood from <a href="http://www.geometrie.tugraz.at/wallner/sgp06-electronic.pdf" target="_blank">Yang et al.</a> are used.

The algorithm automatically detects a suitable radius for each vertex corresponding to its scale and deals with local noise to generate a robust curvature estimate.

Note: This is research code and not really optimized for runtime / cross-platform usage.

## Building

    $ git clone https://github.com/seepa/curvest
    $ cd curvest
    $ mkdir build && cd $_
    $ cmake .. && make

#### Requirements

* Linux, gcc (never tested on windows)
* cmake >= 3.0
* gnuplot (optional)

The code depends on [MVE](https://github.com/simonfuhrmann/mve), which will be automatically downloaded to the build directory.

## Example

Using the example mesh `example/lionhead_small.ply`,
![Lionhead shaded](https://raw.github.com/seepa/curvest/master/example/lion_shaded.png)

computing its mean curvature field is as simple as issuing the following command (assuming you are in the root of this repository):

    $ build/curvestimate example/lionhead_small.ply lionhead_result


The resulting curvature field is saved as *vertex values* in the output .ply file and
can be colorized using the **mesh_colorize** tool (see further down). The result looks like this:
![Lionhead shaded](https://raw.github.com/seepa/curvest/master/example/lion_curvature.png)



## Included Tools

### Curvature Estimation

This is the main application to compute a multi-scale curvature field given a triangle mesh.

    $ ./curvestimate IN-PLY OUT

For a complete list of arguments, run:

    $ ./curvestimate

    Multi-Scale Mean Curvature Estimation

    Usage: ./curvestimate [ OPTIONS ] IN-PLY OUT
    Available options:
      --subdiv=ARG          Number of subdivision iterations for refining the sphere [2]
      --factor-initial=ARG  Initial radius factor [1]
      --factor-inc=ARG      Radius factor increment [1.3]
      --factor-max=ARG      Max radius factor [10]
      --radius-smooth=ARG   Initial radius smoothing iterations [5]
      -e, --edge-smooth=ARG  Edge smoothing factor (0: no smooth, 1.0: full smooth) [0.2]
      -p, --planar-thres=ARG  Planar threshold [0.2]
      -d, --debug-verts=ARG  (Debug) Comma separated list of vertices to debug
      --save-graphs         (Debug) Save a graph for each vertex
      --save-volume-all     (Debug) Save volume at each radius
      --save-patch-all      (Debug) Save patch at each radius
      --save-patch-final    (Debug) Save patch at final radius
      --fixed-radius=ARG    Compute curvature field with a fixed ball radius.
      --use-confidences     Use mesh confidences as final radii.


**Useful arguments**:

* ***-p*** : increase this, if your model is very noisy and you want more smoothing
* ***-e*** : increase this, if you want to smooth edges more
* ***--factor-initial*** : increase this, if you want more global smoothing (the start radius at each vertex will be bigger)


### Mesh Colorize

A small helper application which uses [libcolormap](https://github.com/seepa/colormap), an implementation of diverging colormaps, to colorize a mesh from its vertex values/confidences.

    Colorize a mesh from its vertex values/confidences using diverging colormaps.

    Usage: ./mesh_colorize [ OPTIONS ] IN-PLY OUT
    Available options:
      -c, --confidences     Use vertex confidences.
      --range=ARG           Range of the colormap (comma-separated)
      --colors=ARG          RGB colors for low and high ends of the colormap (e.g. 1.0,0.0,0.0,0.0,0.0,1.0 for red-blue)
      --colorbar            Save the colorbar as a PNG.

By default vertex values and the red-blue, aka *RdBu*, colormap are used.


### Density Field

Creates a density field from the curvature field by remapping the mean curvature values suitable for mesh simplification.

## Possible Applications

Once computed, the curvature field can be used for, e.g.:

* **Mesh simplification** or
* **Mesh smoothing**

In both cases, the curvature field may guide the algorithm in a way which preserves features of the model while smoothing (noisy) planar regions.

## About

This project resulted from my Bachelor's Thesis *Multi-Scale Curvature Field of Triangle Meshes* which I handed in and defended in January 2016 at TU Darmstadt.

TODO: add a link to the thesis/paper

## License

The source code is licensed under the BSD 3-Clause, see `LICENSE.txt`.
