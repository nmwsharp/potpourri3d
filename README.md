# potpourri3d

A Python library of various algorithms and utilities for 3D triangle meshes, polygon meshes, and point clouds. Managed by [Nicholas Sharp](https://nmwsharp.com), with new tools added lazily as needed. Currently, mainly bindings to C++ tools from [geometry-central](http://geometry-central.net/).

`pip install potpourri3d`

The blend includes:
- Mesh and point cloud reading/writing to a few file formats
- Use **heat methods** to compute unsigned and signed distances, parallel transport, logarithmic maps, and more
- Computing geodesic polylines along surface via edge flips
- More!

## Installation

Potpourri3d is on the pypi package index with precompiled binaries for most configuations. Get it like:

`pip install potpourri3d`

If none of the precompiled binaries match your system, `pip` will attempt to compile the library from scratch. This requires `cmake` and a workng C++ compiler toolchain.

**Note**: Some bound functions invoke sparse linear solvers internally. The precompiled binaries use Eigen's solvers; using Suitesparse's solvers instead may significantly improve performance & robustness. To get them, locally compile the package on a machine with Suitesparse installed using the command below ([relevant docs](http://geometry-central.net/build/dependencies/#suitesparse)).

```
python -m pip install potpourri3d --no-binary potpourri3d
```

## Documentation

- [Input / Output](#input--output)
- [Mesh basic utilities](#mesh-basic-utilities)
- [Mesh Distance](#mesh-distance)
- [Mesh Signed Distance](#mesh-signed-distance)
- [Mesh Fast Marching Distance](#mesh-fast-marching-distance)
- [Mesh Vector Heat](#mesh-vector-heat)
- [Mesh Geodesic Paths](#mesh-geodesic-paths)
- [Mesh Geodesic Tracing](#mesh-geodesic-tracing)
- [Polygon Mesh Distance & Transport](#polygon-mesh-distance--transport)
- [Point Cloud Distance & Vector Heat](#point-cloud-distance--vector-heat)
- [Other Point Cloud Routines](#other-point-cloud-routines)

### Input / Output

Read/write meshes and point clouds from some common formats.

- `read_mesh(filename)` Reads a mesh from file. Returns numpy matrices `V, F`, a Nx3 real numpy array of vertices and a Mx3 integer numpy array of 0-based face indices (or Mx4 for a quad mesh, etc).
- `read_polygon_mesh(filename)` Reads a mesh from file. Returns numpy matrices `V, F`, where `V` is a Nx3 real numpy array of vertices, and a `polygons` is a nested list of integers; each sub-list represents a polygon face with 0-based face indices.
  - `filename` the path to read the file from. Currently supports the same file types as [geometry-central](http://geometry-central.net/surface/utilities/io/#supported-file-types). The file type is inferred automatically from the path extension.

- `write_mesh(V, F, filename, UV_coords=None, UV_type=None)` Write a mesh to file, optionally with UV coords.
  - `V` a Nx3 real numpy array of vertices 
  - `F` a Mx3 integer numpy array of faces, with 0-based vertex indices  (or Mx4 for a quad mesh, etc).
  - `filename` the path to write the file to. Currently supports the same file types as [geometry-central](http://geometry-central.net/surface/utilities/io/#supported-file-types). The file type is inferred automatically from the path extension.
  - `UV_coords` (optional) a Ux2 numpy array of UV coords, interpreted based on UV_type. *Warning:* this function does not currently preserve shared UV indices when writing, each written coordinate is independent
  - `UV_type` (optional) string, one of `'per-vertex'`, `'per-face'`, or `'per-corner'`. The size of `U` should be `N`, `M`, or `M*3/4`, respectively

- `read_point_cloud(filename)` Reads a point cloud from file. Returns numpy matrix `V`, a Nx3 real numpy array of vertices.  Really, this just reads a mesh file and ignores the face entries.
  - `filename` the path to read the file from. Currently supports the same file types as [geometry-central](http://geometry-central.net/surface/utilities/io/#supported-file-types)'s mesh reader. The file type is inferred automatically from the path extension.

- `write_point_cloud(V, filename)` Write a mesh to file. Really, this just writes a mesh file with no face entries.
  - `V` a Nx3 real numpy array of vertices 
  - `filename` the path to write the file to. Currently supports the same file types as [geometry-central](http://geometry-central.net/surface/utilities/io/#supported-file-types)'s mesh writer. The file type is inferred automatically from the path extension.

### Mesh basic utilities

- `face_areas(V, F)` computes a length-F real numpy array of face areas for a triangular mesh
- `vertex_areas(V, F)` computes a length-V real numpy array of vertex areas for a triangular mesh (equal to 1/3 the sum of the incident face areas)
- `cotan_laplacian(V, F, denom_eps=0.)` computes the cotan-Laplace matrix as a VxV real sparse csr scipy matrix. Optionally, set `denom_eps` to a small value like `1e-6` to get some additional stability in the presence of degenerate faces.
- `edges(V, F)` returns the Ex2 integer-valued matrix representing the edges of the given surface mesh, as constructed internally. The i-th row gives the indices of the i-th edge's two endpoint vertices.
- _Barycentric points_ are used as input and output to some algorithms below, specified as 2-tuples of the form `(element_index, barycentric_coordinates)`:
  - Vertices are specified as `(vertex_index, )`
  - Edges are specified as `(edge_index, [t])` where t ∈ [0,1] is the parameter along the edge
  - Faces are specified as `(face_index, [tA, tB])` where `tA`, `tB` (and optionally, `tC`) are barycentric coordinates in the face. If `tC` is not specified, then `tC` is inferred to be `1 - tA - tB`.
- `MarchingTrianglesSolver(V, F)` construct an instance of a solver class for contouring scalar functions on triangle meshes using the marching triangles algorithm. 
  - `MarchingTrianglesSolver.marching_triangles(u, isoval=0.)` takes as input a vector `u` representing a scalar function defined on mesh vertices, and an isovalue; returns a list of lists of barycentric points, where each sublist represents a single connected curve component.
  - `marching_triangles(V, F, u, isoval=0.)` is similar to the above, but one-off instead of stateful. Returns a list of lists of barycentric points.

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

The heat method works by solving a sequence of linear PDEs on the surface of your shape. On extremely coarse meshes, it may yield inaccurate results, if you observe this, consider using a finer mesh to improve accuracy. (TODO: do this internally with intrinsic Delaunay refinement.)

- `MeshHeatMethodDistanceSolver(V, F, t_coef=1., use_robust=True)` construct an instance of the solver class.
  - `V` a Nx3 real numpy array of vertices 
  - `F` a Mx3 integer numpy array of faces, with 0-based vertex indices (triangle meshes only, but need not be manifold).
  - `t_coef` set the time used for short-time heat flow. Generally don't change this. If necessary, larger values may make the solution more stable at the cost of smoothing it out.
  - `use_robust` use intrinsic triangulations for increased robustness. Generaly leave this enabled.  
- `MeshHeatMethodDistanceSolver.compute_distance(v_ind)` compute distance from a single vertex, given by zero-based index. Returns an array of distances.
- `MeshHeatMethodDistanceSolver.compute_distance_multisource(v_ind_list)` compute distance from the nearest of a collection of vertices, given by a list of zero-based indices. Returns an array of distances.
- `compute_distance(V, F, v_ind)` Similar to above, but one-off instead of stateful. Returns an array of distances.
- `compute_distance_multisource(V, F, v_ind_list)` Similar to above, but one-off instead of stateful. Returns an array of distances.

### Mesh Signed Distance

Use the [signed heat method](https://nzfeng.github.io/research/SignedHeatMethod/index.html) to compute signed distance on meshes, robust to holes and noise. Repeated solves are fast after initial setup.

```python
import potpourri3d as pp3d

V, F = # your mesh
solver = pp3d.MeshSignedHeatSolver(V, F)

# Specify a curve as a sequence of barycentric points
curves = [
           [
             (61, [0.3, 0.3]), # face
             (7, []), # vertex
             (16, [0.3, 0.3, 0.4]), # face
             (11, [0.4]), # edge
             (71, []), # vertex
             (20, [0.3, 0.3, 0.4]), # face
             (13, []), # vertex
             (58, []) # vertex
             ]
         ]

# Compute a distance field combining signed distance to curve sources, and unsigned distance to point sources.
dist = solver.compute_distance(curves, [], points) 
```

- `MeshSignedHeatSolver.compute_distance(curves, is_signed, points, preserve_source_normals=False, level_set_constraint="ZeroSet", soft_level_set_weight=-1)`
  - `curves` a list of lists of source points; each point is specified via barycentric coordinates.
  - `is_signed` a list of bools, one for each curve in `curves`, indicating whether one should compute signed distance (`True`) or unsigned distance (`False`) to a curve. All `True` by default.
  - `points` a list of source vertex indices
  - `preserve_source_normals` whether to additionally constrain the normals of the curve. Generally not necessary.
  - `level_set_constraint` whether to apply level set constraints, with options "ZeroSet", "None", "Multiple". Generally set to "ZeroSet" (set by default).
  - `soft_level_set_weight` float; if positive, gives the weight with which the given level set constraint is "softly" enforced (negative by default). Generally not necessary.

### Mesh Fast Marching Distance

```python
import potpourri3d as pp3d

V, F = # your mesh
solver = pp3d.MeshFastMarchingDistanceSolver(V, F)

# Specify each curve as a sequence of barycentric points
curves = [
           [
             (61, [0.3, 0.3]), # face
             (7, []), # vertex
             (16, [0.3, 0.3, 0.4]), # face
             (11, [0.4]), # edge
             (71, []), # vertex
             (20, [0.3, 0.3, 0.4]), # face
             (13, []), # vertex
             (58, []) # vertex
             ]
         ]

# Compute a signed distance field to a set of closed curves.
signed_dist = solver.compute_distance(curves, sign=True) 

# Compute unsigned to a set of points.
points = [
          [
            (71, []), # vertex
            (18, [0.5]) # edge
          ]
         ]
unsigned_dist = solver.compute_distance(points, sign=False) 
```

- `MeshFastMarchingDistanceSolver.compute_distance(curves, distances=[], sign=False)`
  - `curves` a list of lists of source points; each point is specified via barycentric coordinates.
  - `distances` a list of lists of initial distances. Default initial distances are 0.
  - `sign` if False, compute unsigned distance; if True, compute signed distance. (When initial distances are not 0, "signed" means that the gradient of distance is continuous across the source curves.)

### Mesh Vector Heat

Use the [vector heat method](https://nmwsharp.com/research/vector-heat-method/) and [affine heat method](https://www.yousufsoliman.com/projects/the-affine-heat-method.html) to compute various interpolation & vector-based quantities on meshes. Repeated solves are fast after initial setup.

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
```

- `MeshVectorHeatSolver(V, F, t_coef=1., useIntrinsicDelaunay=True)` construct an instance of the solver class.
  - `V` a Nx3 real numpy array of vertices 
  - `F` a Mx3 integer numpy array of faces, with 0-based vertex indices (triangle meshes only, should be manifold).
  - `t_coef` set the time used for short-time heat flow. Generally don't change this. If necessary, larger values may make the solution more stable at the cost of smoothing it out.
  - `useIntrinsicDelaunay` if true, an [intrinsic triangulation](https://geometry-central.net/surface/intrinsic_triangulations/basics/) is used internally to improve robustness
- `MeshVectorHeatSolver.extend_scalar(v_inds, values)` nearest-geodesic-neighbor interpolate values defined at vertices. Vertices will take the value from the closest source vertex (plus some slight smoothing)
  - `v_inds` a list of source vertices
  - `values` a list of scalar values, one for each source vertex
- `MeshVectorHeatSolver.get_tangent_frames()` get the coordinate frames used to define tangent data at each vertex. Returned as a tuple of basis-X, basis-Y, and normal axes, each as an Nx3 array. May be necessary for change-of-basis into or out of tangent vector convention.
- `MeshVectorHeatSolver.get_connection_laplacian()` get the _connection Laplacian_ used internally in the vector heat method, as a VxV sparse matrix.
- `MeshVectorHeatSolver.transport_tangent_vector(v_ind, vector)` parallel transports a single vector across a surface
  - `v_ind` index of the source vertex
  - `vector` a 2D tangent vector to transport
- `MeshVectorHeatSolver.transport_tangent_vectors(v_inds, vectors)` parallel transports a collection of vectors across a surface, such that each vertex takes the vector from its nearest-geodesic-neighbor.
  - `v_inds` a list of source vertices
  - `vectors` a list of 2D tangent vectors, one for each source vertex
- `MeshVectorHeatSolver.compute_log_map(v_ind, strategy='AffineLocal')` compute the logarithmic map centered at the given source vertex
  - `v_ind` index of the source vertex
  - `strategy` one of `'VectorHeat'`,`'AffineLocal'`, `'AffineAdaptive'`, see [here for an explanation](https://geometry-central.net/surface/algorithms/vector_heat_method/#logarithmic-map)


### Mesh Geodesic Paths

Use [edge flips to compute geodesic paths](https://nmwsharp.com/research/flip-geodesics/) on surfaces. These methods take an initial path, loop, or start & end points along the surface, and straighten the path out to be geodesic.

This approach is mainly useful when you want the path itself, rather than the distance. These routines use an iterative strategy which is quite fast, but note that it is not guaranteed to generate a globally-shortest geodesic (they sometimes find some other very short geodesic instead if straightening falls into different local minimum).

<img src="https://github.com/nmwsharp/potpourri3d/blob/master/media/elephant_geodesic.jpg" height="400">

```python
import potpourri3d as pp3d

V, F = # your mesh
path_solver = pp3d.EdgeFlipGeodesicSolver(V,F) # shares precomputation for repeated solves
path_pts = path_solver.find_geodesic_path(v_start=14, v_end=22)
# path_pts is a Vx3 numpy array of points forming the path
```

- `EdgeFlipGeodesicSolver(V, F)` construct an instance of the solver class.
  - `V` a Nx3 real numpy array of vertices 
  - `F` a Mx3 integer numpy array of faces, with 0-based vertex indices (must form a manifold, oriented triangle mesh).
- `EdgeFlipGeodesicSolver.find_geodesic_path(v_start, v_end, max_iterations=None, max_relative_length_decrease=None)` compute a geodesic from `v_start` to `v_end`. Output is an `Nx3` numpy array of positions which define the path as a polyline along the surface.
- `EdgeFlipGeodesicSolver.find_geodesic_path_poly(v_list, max_iterations=None, max_relative_length_decrease=None)` like `find_geodesic_path()`, but takes as input a list of vertices `[v_start, v_a, v_b, ..., v_end]`, which is shorted to find a path from `v_start` to `v_end`. Useful for finding geodesics which are not shortest paths. The input vertices do not need to be connected; the routine internally constructs a piecwise-Dijkstra path between them. However, that path must not cross itself.
- `EdgeFlipGeodesicSolver.find_geodesic_loop(v_list, max_iterations=None, max_relative_length_decrease=None)` like `find_geodesic_path_poly()`, but connects the first to last point to find a closed geodesic loop.

In the functions above, the optional argument `max_iterations` is an integer, giving the the maximum number of shortening iterations to perform (default: no limit). The optional argument `max_relative_length_decrease` is a float limiting the maximum decrease in length for the path, e.g. `0.5` would mean the resulting path is at least `0.5 * L` length, where `L` is the initial length.

### Mesh Geodesic Tracing

Given an initial point and direction/length, these routines trace out a geodesic path along the surface of the mesh and return it as a polyline.

```python
import potpourri3d as pp3d

V, F = # your mesh
tracer = pp3d.GeodesicTracer(V,F) # shares precomputation for repeated traces

trace_pts = tracer.trace_geodesic_from_vertex(22, np.array((0.3, 0.5, 0.4)))
# trace_pts is a Vx3 numpy array of points forming the path
```

- `GeodesicTracer(V, F)` construct an instance of the tracer class.
  - `V` a Nx3 real numpy array of vertices
  - `F` a Mx3 integer numpy array of faces, with 0-based vertex indices (must form a manifold, oriented triangle mesh).
- `GeodesicTracer.trace_geodesic_from_vertex(start_vert, direction_xyz, max_iterations=None)` trace a geodesic from `start_vert`. `direction_xyz` is a length-3 vector giving the direction to walk trace in 3D xyz coordinates, it will be projected onto the tangent space of the vertex. The magnitude of `direction_xyz` determines the distance walked. Output is an `Nx3` numpy array of positions which define the path as a polyline along the surface.
- `GeodesicTracer.trace_geodesic_from_face(start_face, bary_coords, direction_xyz, max_iterations=None)` similar to above, but from a point in a face. `bary_coords` is a length-3 vector of barycentric coordinates giving the location within the face to start from.

Set `max_iterations` to terminate early after tracing the path through some number of faces/edges (default: no limit).

### Polygon Mesh Distance & Transport

Use the [heat method for unsigned geodesic distance](https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/), the [signed heat method](https://nzfeng.github.io/research/SignedHeatMethod/index.html) to compute signed distance, and the [vector heat method](https://nmwsharp.com/research/vector-heat-method/) to compute various interpolation & vector-based quantities on general polygon meshes (including mixed-degree meshes, such as tri-quad meshes). Repeated solves are fast after initial setup.

```python
import potpourri3d as pp3d

V, polygons = # your polygon mesh
solver = pp3d.PolygonMeshHeatSolver(V, F)

# Compute unsigned geodesic distance to vertices 12 and 17
dist = solver.compute_distance([12, 17])

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

# Compute signed distance to the oriented curve(s) denoted by a vertex sequence.
curves = [
           [9, 10, 12, 13, 51, 48], 
           [79, 93, 12, 30, 78, 18, 92], 
           [90, 84, 19, 91, 82, 81, 83]
         ]
signed_dist = solver.compute_signed_distance(curves)
```

- `PolygonMeshHeatSolver(V, polygons, t_coef=1.)` construct an instance of the solver class.
  - `V` a Nx3 real numpy array of vertices 
  - `polygons` a list of lists; each sub-list represents a polygon face with 0-based face indices (integers).
  - `t_coef` set the time used for short-time heat flow. Generally don't change this. If necessary, larger values may make the solution more stable at the cost of smoothing it out.
- `PolygonMeshHeatSolver.extend_scalar(v_inds, values)` nearest-geodesic-neighbor interpolate values defined at vertices. Vertices will take the value from the closest source vertex (plus some slight smoothing)
  - `v_inds` a list of source vertices
  - `values` a list of scalar values, one for each source vertex
- `PolygonMeshHeatSolver.get_tangent_frames()` get the coordinate frames used to define tangent data at each vertex. Returned as a tuple of basis-X, basis-Y, and normal axes, each as an Nx3 array. May be necessary for change-of-basis into or out of tangent vector convention.
- `PolygonMeshHeatSolver.transport_tangent_vectors(v_inds, vectors)` parallel transports a collection of vectors across a surface, such that each vertex takes the vector from its nearest-geodesic-neighbor.
  - `v_inds` a list of source vertices
  - `vectors` a list of 2D tangent vectors, one for each source vertex
- `PolygonMeshHeatSolver.compute_distance(v_inds)` 
  - `v_inds` a list of source vertices
- `PolygonMeshHeatSolver.compute_signed_distance(curves, level_set_constraint="ZeroSet")`
  - `curves` a list of lists of source vertices
  - `level_set_constraint` whether to apply level set constraints, with options "ZeroSet", "None", "Multiple". Generally set to "ZeroSet" (set by default).

 
### Point Cloud Distance & Vector Heat

Use the [heat method for unsigned geodesic distance](https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/), the [signed heat method](https://nzfeng.github.io/research/SignedHeatMethod/index.html) to compute signed distance, and the [vector heat method](https://nmwsharp.com/research/vector-heat-method/) to compute various interpolation & vector-based quantities on point clouds. Repeated solves are fast after initial setup.

![point cloud vector heat examples](https://github.com/nmwsharp/potpourri3d/blob/master/media/point_heat_solvers.jpg)

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

# Signed distance to the oriented curve(s) denoted by a point sequence.
curves = [
           [9, 10, 12, 13, 51, 48], 
           [79, 93, 12, 30, 78, 18, 92], 
           [90, 84, 19, 91, 82, 81, 83]
         ]
signed_dist = solver.compute_signed_distance(curves, basisN)
```

- `PointCloudHeatSolver(P, t_coef=1.)` construct an instance of the solver class.
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
- `PointCloudHeatSolver.compute_signed_distance(curves, cloud_normals, preserve_source_normals=False, level_set_constraint="ZeroSet", soft_level_set_weight=-1)`
  - `curves` a list of lists of source point indices
  - `cloud_normals` a list of 3D normal vectors, one for each point in the point cloud
  - `preserve_source_normals` whether to additionally constrain the normals of the curve. Generally not necessary.
  - `level_set_constraint` whether to apply level set constraints, with options "ZeroSet", "None", "Multiple". Generally set to "ZeroSet" (set by default).
  - `soft_level_set_weight` float; if positive, gives the weight with which the given level set constraint is "softly" enforced (negative by default). Generally not necessary.


### Other Point Cloud Routines

#### Local Triangulation

Construct a _local triangulation_ of a point cloud, a surface-like set of triangles amongst the points in the cloud. This is _not_ a nice connected/watertight mesh, instead it is a crazy soup, which is a union of sets of triangles computed independently around each point. These triangles _are_ suitable for running many geometric algorithms on, such as approximating surface properties of the point cloud, evaluating physical and geometric energies, or building Laplace matrices. See "A Laplacian for Nonmanifold Triangle Meshes", Sharp & Crane 2020, Sec 5.7 for details.

- `PointCloudLocalTriangulation(P, with_degeneracy_heuristic=True)`
  - `PointCloudLocalTriangulation.get_local_triangulation()` returns a `[V,M,3]` integer numpy array, holding indices of vertices which form the triangulation. Each `[i,:,:]` holds the local triangles about vertex `i`. `M` is the max number of neighbors in any local triangulation. For vertices with fewer neighbors, the trailing rows hold `-1`.  
