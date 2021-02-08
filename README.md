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

- `read_mesh(filename)` Reads a mesh from file. Returns numpy matrices `V, F`, a Nx3 real numpy array of vertices and a Mx3 integer numpy array of 0-based face indices (or Mx4 for a quad mesh, etc).
  - `filename` the path to read the file from. Currently supports the same file types as [geometry-central](http://geometry-central.net/surface/utilities/io/#supported-file-types). The file type is inferred automatically from the path extension.

- `write_mesh(V, F, filename)` Write a mesh from file. Returns numpy matrices `V, F`, a Vx3 real array of vertices and a Fx3 integer array of 0-based face indices (or Fx4 for a quad mesh, etc).
  - `V` a Nx3 real numpy array of vertices 
  - `F` a Mx3 integer numpy array of 0-based face indices (or Mx4 for a quad mesh, etc).
  - `filename` the path to write the file to. Currently supports the same file types as [geometry-central](http://geometry-central.net/surface/utilities/io/#supported-file-types). The file type is inferred automatically from the path extension.


