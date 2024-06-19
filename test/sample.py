import os, sys

import polyscope as ps
import numpy as np
import scipy 
import scipy.sparse
import scipy.sparse.linalg

# Path to where the bindings live
sys.path.append(os.path.join(os.path.dirname(__file__), "../build/"))
sys.path.append(os.path.join(os.path.dirname(__file__), "../src/"))

import potpourri3d as pp3d

ps.init()

# Read input

## = Mesh test
V, F = pp3d.read_mesh("bunny_small.ply")
ps_mesh = ps.register_surface_mesh("mesh", V, F)

# Distance
dists = pp3d.compute_distance(V, F, 4)
ps_mesh.add_scalar_quantity("dist", dists)

# Vector heat
solver = pp3d.MeshVectorHeatSolver(V, F)

# Vector heat (extend scalar)
ext = solver.extend_scalar([1, 22], [0., 6.]) 
ps_mesh.add_scalar_quantity("ext", ext)

# Vector heat (tangent frames)
basisX, basisY, basisN = solver.get_tangent_frames()
ps_mesh.add_vector_quantity("basisX", basisX)
ps_mesh.add_vector_quantity("basisY", basisY)
ps_mesh.add_vector_quantity("basisN", basisN)

# Vector heat (transport vector)
ext = solver.transport_tangent_vector(1, [6., 6.])
ext3D = ext[:,0,np.newaxis] * basisX +  ext[:,1,np.newaxis] * basisY
ps_mesh.add_vector_quantity("transport vec", ext3D)

ext = solver.transport_tangent_vectors([1, 22], [[6., 6.], [3., 4.]])
ext3D = ext[:,0,np.newaxis] * basisX +  ext[:,1,np.newaxis] * basisY
ps_mesh.add_vector_quantity("transport vec2", ext3D)

# Vector heat (log map)
logmap = solver.compute_log_map(1)
ps_mesh.add_parameterization_quantity("logmap", logmap)

# Flip geodesics
path_solver = pp3d.EdgeFlipGeodesicSolver(V,F)
for k in range(50):
    for i in range(5):
        path_pts = path_solver.find_geodesic_path(v_start=1, v_end=22+i)
        ps.register_curve_network("flip path " + str(i), path_pts, edges='line')
        
path_pts = path_solver.find_geodesic_path_poly([1173, 148, 870, 898])
ps.register_curve_network("flip path poly", path_pts, edges='line')

loop_pts = path_solver.find_geodesic_loop([1173, 148, 870, 898])
ps.register_curve_network("flip loop", loop_pts, edges='loop')

loop_pts = path_solver.find_geodesic_loop([307, 757, 190], max_relative_length_decrease=.9) # this one otherwise contracts to a point
ps.register_curve_network("flip loop", loop_pts, edges='loop')


## = Point cloud test
P = V
ps_cloud = ps.register_point_cloud("cloud", P)

# == heat solver
solver = pp3d.PointCloudHeatSolver(P)

# distance
dists = solver.compute_distance(4)
dists2 = solver.compute_distance_multisource([4, 13, 784])
ps_cloud.add_scalar_quantity("dist", dists)
ps_cloud.add_scalar_quantity("dist2", dists2)

# scalar extension
ext = solver.extend_scalar([1, 22], [0., 6.])
ps_cloud.add_scalar_quantity("ext", ext)

# Vector heat (tangent frames)
basisX, basisY, basisN = solver.get_tangent_frames()
ps_cloud.add_vector_quantity("basisX", basisX)
ps_cloud.add_vector_quantity("basisY", basisY)
ps_cloud.add_vector_quantity("basisN", basisN)

# Vector heat (transport vector)
ext = solver.transport_tangent_vector(1, [6., 6.])
ext3D = ext[:,0,np.newaxis] * basisX +  ext[:,1,np.newaxis] * basisY
ps_cloud.add_vector_quantity("transport vec", ext3D)

ext = solver.transport_tangent_vectors([1, 22], [[6., 6.], [3., 4.]])
ext3D = ext[:,0,np.newaxis] * basisX +  ext[:,1,np.newaxis] * basisY
ps_cloud.add_vector_quantity("transport vec2", ext3D)

# Vector heat (log map)
logmap = solver.compute_log_map(1)
ps_cloud.add_scalar_quantity("logmapX", logmap[:,0])
ps_cloud.add_scalar_quantity("logmapY", logmap[:,1])

# Areas
vert_area = pp3d.vertex_areas(V,F)
ps_mesh.add_scalar_quantity("vert area", vert_area)
face_area = pp3d.face_areas(V,F)
ps_mesh.add_scalar_quantity("face area", face_area, defined_on='faces')


# Laplacian
L = pp3d.cotan_laplacian(V,F,denom_eps=1e-6)
M = scipy.sparse.diags(vert_area)
k_eig = 6
evals, evecs = scipy.sparse.linalg.eigsh(L, k_eig, M, sigma=1e-8)
for i in range(k_eig):
    ps_mesh.add_scalar_quantity("evec " + str(i), evecs[:,i])

ps.show()
