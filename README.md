# potpourri3d

A Python library of various algorithms and utilities for 3D triangle meshes and point clouds. Currently, mainly bindings to C++ tools from [geometry-central](http://geometry-central.net/).

`pip install potpourri3d`

The blend includes:
- Mesh and point cloud reading/writing to a few file formats
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
solver = pp3d.MeshHeatMethodDistanceSolver(V,F)
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

### Mesh Vector Heat

Use the [vector heat method](https://nmwsharp.com/research/vector-heat-method/) to compute various interpolation & vector-based quantities on meshes. Repeated solves are fast after initial setup.

```python
import potpourri3d as pp3d

# = Stateful solves
V, F = # a Nx3 numpy array of points and Mx3 array of triangle face indices
solver = pp3d.MeshVectorHeatSolver(V,F)

# Extend the value `0.` from vertex 12 and `1.` from vertex 17. Any vertex 
# geodesically closer to 12. will take the value 0., and vice versa 
# (plus some slight smoothing)
ext = solver.extend_scalar([12, 17], [0.,1.])

# Get the tangent frames which are used by the solver to define tangent data
# at each vertex
basisX, basisY, basisN = solver.get_tangent_frames()

# Parallel transport a vector along the surface
# (and map it to a vector in 3D)
sourceV = 22
ext = solver.transport_tangent_vector(sourceV, [6., 6.])
ext3D = ext[:,0,np.newaxis] * basisX +  ext[:,1,np.newaxis] * basisY

# Compute the logarithmic map
logmap = solver.compute_log_map(sourceV)
ps_mesh.add_parameterization_quantity("logmap", logmap)
```

- `MeshVectorHeatSolver(self, V, F, t_coef=1.)` construct an instance of the solver class.
  - `V` a Nx3 real numpy array of vertices 
  - `F` a Mx3 integer numpy array of faces, with 0-based vertex indices (triangle meshes only, should be manifold).
  - `t_coef` set the time used for short-time heat flow. Generally don't change this. If necessary, larger values may make the solution more stable at the cost of smoothing it out.
- `MeshVectorHeatSolver.extend_scalar(v_inds, values)` nearest-geodesic-neighbor interpolate values defined at vertices. Vertices will take the value from the closest source vertex (plus some slight smoothing)
  - `v_inds` a list of source vertices
  - `values` a list of scalar values, one for each source vertex
- `MeshVectorHeatSolver.get_tangent_frames()` get the coordinate frames used to define tangent data at each vertex. Returned as a tuple of basis-X, basis-Y, and normal axes, each as an Nx3 array. May be necessary for change-of-basis into or out of tangent vector convention.
- `MeshVectorHeatSolver.transport_tangent_vector(v_ind, vector)` parallel transports a single vector across a surface
  - `v_ind` index of the source vertex
  - `vector` a 2D tangent vector to transport
- `MeshVectorHeatSolver.transport_tangent_vectors(v_inds, vectors)` parallel transports a collection of vectors across a surface, such that each vertex takes the vector from its nearest-geodesic-neighbor.
  - `v_inds` a list of source vertices
  - `vectors` a list of 2D tangent vectors, one for each source vertex
- `MeshVectorHeatSolver.compute_log_map(v_ind)` compute the logarithmic map centered at the given source vertex
  - `v_ind` index of the source vertex

### Point Cloud Distance & Vector Heat

Use the [heat method for geodesic distance](https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/) and [vector heat method](https://nmwsharp.com/research/vector-heat-method/) to compute various interpolation & vector-based quantities on point clouds. Repeated solves are fast after initial setup.


```python
import potpourri3d as pp3d

# = Stateful solves
P = # a Nx3 numpy array of points
solver = pp3d.PointCloudHeatSolver(P)

# Compute the geodesic distance to point 4
dists = solver.compute_distance(4)

# Extend the value `0.` from point 12 and `1.` from point 17. Any point 
# geodesically closer to 12. will take the value 0., and vice versa 
# (plus some slight smoothing)
ext = solver.extend_scalar([12, 17], [0.,1.])

# Get the tangent frames which are used by the solver to define tangent data
# at each point
basisX, basisY, basisN = solver.get_tangent_frames()

# Parallel transport a vector along the surface
# (and map it to a vector in 3D)
sourceP = 22
ext = solver.transport_tangent_vector(sourceP, [6., 6.])
ext3D = ext[:,0,np.newaxis] * basisX +  ext[:,1,np.newaxis] * basisY

# Compute the logarithmic map
logmap = solver.compute_log_map(sourceP)
```

- `PointCloudHeatSolver(self, P, t_coef=1.)` construct an instance of the solver class.
  - `P` a Nx3 real numpy array of points
  - `t_coef` set the time used for short-time heat flow. Generally don't change this. If necessary, larger values may make the solution more stable at the cost of smoothing it out.
- `PointCloudHeatSolver.extend_scalar(p_inds, values)` nearest-geodesic-neighbor interpolate values defined at points. Points will take the value from the closest source point (plus some slight smoothing)
  - `v_inds` a list of source points
  - `values` a list of scalar values, one for each source points
- `PointCloudHeatSolver.get_tangent_frames()` get the coordinate frames used to define tangent data at each point. Returned as a tuple of basis-X, basis-Y, and normal axes, each as an Nx3 array. May be necessary for change-of-basis into or out of tangent vector convention.
- `PointCloudHeatSolver.transport_tangent_vector(p_ind, vector)` parallel transports a single vector across a surface
  - `p_ind` index of the source point
  - `vector` a 2D tangent vector to transport
- `PointCloudHeatSolver.transport_tangent_vectors(p_inds, vectors)` parallel transports a collection of vectors across a surface, such that each vertex takes the vector from its nearest-geodesic-neighbor.
  - `p_inds` a list of source points
  - `vectors` a list of 2D tangent vectors, one for each source point
- `PointCloudHeatSolver.compute_log_map(p_ind)` compute the logarithmic map centered at the given source point
  - `p_ind` index of the source point


