import os, sys

import polyscope as ps
import numpy as np
# import scipy.sparse.linalg as sla

# Path to where the bindings live
sys.path.append(os.path.join(os.path.dirname(__file__), "../build/"))
sys.path.append(os.path.join(os.path.dirname(__file__), "../src/"))

import potpourri3d as pp3d

ps.init()

# Read input

## = Mesh test
# V, F = pp3d.read_mesh("/Users/nick/mesh/spot.obj")
V, F = pp3d.read_mesh("bunny_small.ply")
dists = pp3d.compute_distance(V, F, 4)

solver = pp3d.MeshVectorHeatSolver(V, F)
ext = solver.extend_scalar([1, 22], [0., 6.])

ps_mesh = ps.register_surface_mesh("mesh", V, F)
ps_mesh.add_scalar_quantity("dist", dists)
ps_mesh.add_scalar_quantity("ext", ext)


## = Point cloud test
P = V

solver = pp3d.PointCloudHeatSolver(P)
dists = solver.compute_distance(4)
dists2 = solver.compute_distance_multisource([4, 13, 784])
ext = solver.extend_scalar([1, 22], [0., 6.])

ps_cloud = ps.register_point_cloud("cloud", P)
ps_cloud.add_scalar_quantity("dist", dists)
ps_cloud.add_scalar_quantity("dist2", dists2)
ps_cloud.add_scalar_quantity("ext", ext)

ps.show()
