# potpourri3d

A Python library of various algorithms and utilities for 3D triangle meshes and point clouds. Currently, mainly bindings to C++ tools from [geometry-central](http://geometry-central.net/).

`pip install potpourri3d`

The blend includes:
- Mesh and point cloud reading/writing to a few file formats
- Sample a point cloud from a mesh
- Use **heat methods** to compute distance, parallel transport, logarithmic maps, and more

## Installation

Potpourri3d is on the pypi package index with precompiled binaries for most configuations. Get it like:

`pip install potpourri3d`

If none of the precompiled binaries match your system, `pip` will attempt to compile the library from scratch. This requires `cmake` and a workng C++ compiler toolchain.

**Note**: Some bound functions invoke sparse linear solvers internally. The precompiled binaries use Eigen's solvers; using Suitesparse's solvers instead may significantly improve performance & robustness. To get them, locally compile the package on a machine with Suitesparse installed using the command below ([relevant docs](http://geometry-central.net/build/dependencies/#suitesparse)).

```
python -m pip install potpourri3d --no-binary potpourri3d
```

## Documentation

### IO

Read/write meshes and point clouds from some common formats.

- `read_mesh(filename)` Reads a mesh from file. Returns numpy matrices `V, F`, a Nx3 real numpy array of vertices and a Mx3 integer numpy array of 0-based face indices (or Mx4 for a quad mesh, etc).
  - `filename` the path to read the file from. Currently supports the same file types as [geometry-central](http://geometry-central.net/surface/utilities/io/#supported-file-types). The file type is inferred automatically from the path extension.

- `write_mesh(V, F, filename)` Write a mesh from file. Returns numpy matrices `V, F`, a Vx3 real array of vertices and a Fx3 integer array of 0-based face indices (or Fx4 for a quad mesh, etc).
  - `V` a Nx3 real numpy array of vertices 
  - `F` a Mx3 integer numpy array of faces, with 0-based vertex indices  (or Mx4 for a quad mesh, etc).
  - `filename` the path to write the file to. Currently supports the same file types as [geometry-central](http://geometry-central.net/surface/utilities/io/#supported-file-types). The file type is inferred automatically from the path extension.

### Mesh Distance

Use the [heat method for geodesic distance](https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/) to compute geodesic distance on surfaces. Repeated solves are fast after initial setup. Uses [intrinsic triangulations](http://www.cs.cmu.edu/~kmcrane/Projects/NonmanifoldLaplace/NonmanifoldLaplace.pdf) internally for increased robustness.

```python
import potpourri3d as pp3d

# = Stateful solves (much faster if computing distance many times)
solver = pp3d.MeshHeatMethodDistanceSolver(V,F,7)
dist = solver.compute_distance(7)
dist = solver.compute_distance_multisource([1,2,3])  

# = One-off versions
dist = pp3d.compute_distance(V,F,7)
dist = pp3d.compute_distance_multisource(V,F,[1,3,4])
```


- `MeshHeatMethodDistanceSolver(self, V, F, t_coef=1., use_robust=True)` construct an instance of the solver class.
  - `V` a Nx3 real numpy array of vertices 
  - `F` a Mx3 integer numpy array of faces, with 0-based vertex indices (triangle meshes only, but need not be manifold).
  - `t_coef` set the time used for short-time heat flow. Generally don't change this. If necessary, larger values may make the solution more stable at the cost of smoothing it out.
  - `use_robust` use intrinsic triangulations for increased robustness. Generaly leave this enabled.  
- `MeshHeatMethodDistanceSolver.compute_distance(v_ind)` compute distance from a single vertex, given by zero-based index. Returns an array of distances.
- `MeshHeatMethodDistanceSolver.compute_distance_multisource(v_ind_list)` compute distance from the nearest of a collection of vertices, given by a list of zero-based indices. Returns an array of distances.
- `compute_distance(V, F, v_ind)` Similar to above, but one-off instead of stateful. Returns an array of distances.
- `compute_distance_multisource(V, F, v_ind_list)` Similar to above, but one-off instead of stateful. Returns an array of distances.

