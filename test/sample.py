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
ext3D = ext[:,0,np.newaxis] * basisX +	ext[:,1,np.newaxis] * basisY
ps_mesh.add_vector_quantity("transport vec", ext3D)

ext = solver.transport_tangent_vectors([1, 22], [[6., 6.], [3., 4.]])
ext3D = ext[:,0,np.newaxis] * basisX +	ext[:,1,np.newaxis] * basisY
ps_mesh.add_vector_quantity("transport vec2", ext3D)

# Vector heat (log map)
logmap = solver.compute_log_map(1)
ps_mesh.add_parameterization_quantity("logmap", logmap)

# Signed heat
curve_idxs = [809, 808, 830, 783, 843, 181, 696, 781, 302, 937, 793, 26, 1025, 702, 1047, 1088, 487, 981, 1239, 1238, 402, 57, 1189, 161, 471, 529, 560, 316, 766, 833]
broken_curves = [
					 [809, 808, 830, 783, 843, 181, 696, 781, 302, 937, 793],
					 [1088, 487, 981, 1239, 1238, 402, 57, 1189],
					 [529, 560, 316, 766, 833]
				 ]
point_idxs = [1085, 1000]

# Convert to SurfacePoint format
points = [(idx, []) for idx in point_idxs]
curves = [[(idx, []) for idx in curve] for curve in broken_curves]
sp_curves = [
							[(603, [0.3, 0.3, 0.4]), (726, []), (1660, [0.3, 0.3, 0.4]), (718, []), (2058, [0.3, 0.3, 0.4]), (162, []), (583, []), (88, []), (26, []), (1025, []), (600, []), (1394, [0.3, 0.3]), (1101, []), (1517, [0.3, 0.3]), (919, []), (62, [])]
						]

E = pp3d.edges(V, F)
def interpolate_coordinates(curve):
	positions = []
	for point in curve:
		(idx, coords) = point
		if len(coords) < 1:
			positions.append(V[idx,:])
		elif len(coords) == 1:
			pos = (1.-coords[0]) * V[E[idx,0],:] + coords[0] * V[E[idx,1],:]
			positions.append(pos)
		elif len(coords) >= 2:
			pos = coords[0]*V[F[idx,0],:] + coords[1]*V[F[idx,1],:] + (1.-coords[0]-coords[1])*V[F[idx,2],:]
			positions.append(pos)
	return positions

sp_curves_positions = interpolate_coordinates(sp_curves[0])

signed_solver = pp3d.MeshSignedHeatSolver(V, F)
unsigned_dist = signed_solver.compute_distance(curves, [False, False, False], points)
signed_dist = signed_solver.compute_distance(curves, [], [],
	preserve_source_normals=True, level_set_constraint="none", soft_level_set_weight=-1)
sp_dist = signed_solver.compute_distance(sp_curves, [], [],
	preserve_source_normals=False, level_set_constraint="ZeroSet", soft_level_set_weight=-1)
mixed_dist = signed_solver.compute_distance(curves, [True, True, True], points,
	preserve_source_normals=False, level_set_constraint="ZeroSet", soft_level_set_weight=-1)
isocontour = pp3d.marching_triangles(V, F, signed_dist, 0.1)
contour = interpolate_coordinates(isocontour[0])
ps.register_curve_network("broken curve 0", np.array(V[broken_curves[0]]), edges="line")
ps.register_curve_network("broken curve 1", np.array(V[broken_curves[1]]), edges="line")
ps.register_curve_network("broken curve 2", np.array(V[broken_curves[2]]), edges="line")
ps.register_curve_network("curve (barycentric points)", np.array(sp_curves_positions), edges="line")
ps.register_curve_network("isocontour", np.array(contour), edges="line")
ps_mesh.add_distance_quantity("unsigned distance", unsigned_dist)
ps_mesh.add_distance_quantity("signed distance", signed_dist)
ps_mesh.add_distance_quantity("signed distance (barycentric)", sp_dist)
ps_mesh.add_distance_quantity("mixed distance", mixed_dist)

# Fast marching
fmm_solver = pp3d.MeshFastMarchingDistanceSolver(V, F)
fmm_dist = fmm_solver.compute_distance([[(idx, []) for idx in curve_idxs]], sign=True)
ps_mesh.add_distance_quantity("fast marching distance", fmm_dist)

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


# Trace geodesics
tracer = pp3d.GeodesicTracer(V,F)

trace_pts = tracer.trace_geodesic_from_vertex(22, np.array((0.3, 0.5, 0.4)))
ps.register_curve_network("trace vertex geodesic", trace_pts, edges='line')

trace_pts = tracer.trace_geodesic_from_face(31, np.array((0.1, 0.4, 0.5)), np.array((0.3, 0.5, 0.4)))
ps.register_curve_network("trace face geodesic", trace_pts, edges='line')

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

## = Polygon mesh test

V, F = pp3d.read_polygon_mesh("bunny_poly.ply")
ps_polymesh = ps.register_surface_mesh("polygon mesh", V, F)
solver = pp3d.PolygonMeshHeatSolver(V, F)

# Distance
dist = solver.compute_distance([100, 1000])
ps_polymesh.add_distance_quantity("unsigned distance", dist)

# Vector heat (extend scalar)
ext = solver.extend_scalar([100, 1000], [0., 6.]) 
ps_polymesh.add_scalar_quantity("ext", ext)

# Vector heat (transport vector)
basisX, basisY, basisN = solver.get_tangent_frames()
ext = solver.transport_tangent_vectors([100, 1000], [[6., 6.], [3., 4.]])
ext3D = ext[:,0,np.newaxis] * basisX +	ext[:,1,np.newaxis] * basisY
ps_polymesh.add_vector_quantity("transport vec2", ext3D)

# Signed heat 
curves = [
					 [96, 1037, 1238, 1239, 851, 487, 1088, 601, 217, 942, 702, 631, 805, 793, 937, 1237, 302, 781, 185, 923, 908, 843, 909, 910, 889, 809, 833, 766, 317, 267, 296, 345, 289, 484, 465, 985, 1168, 1099, 17, 212]
				 ]
signed_dist = solver.compute_signed_distance(curves)
ps.register_curve_network("polygon mesh curve", np.array(V[curves[0]]), edges="loop")
ps_polymesh.add_distance_quantity("signed distance", signed_dist)

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
ext3D = ext[:,0,np.newaxis] * basisX +	ext[:,1,np.newaxis] * basisY
ps_cloud.add_vector_quantity("transport vec", ext3D)

ext = solver.transport_tangent_vectors([1, 22], [[6., 6.], [3., 4.]])
ext3D = ext[:,0,np.newaxis] * basisX +	ext[:,1,np.newaxis] * basisY
ps_cloud.add_vector_quantity("transport vec2", ext3D)

# Vector heat (log map)
logmap = solver.compute_log_map(1)
ps_cloud.add_scalar_quantity("logmapX", logmap[:,0])
ps_cloud.add_scalar_quantity("logmapY", logmap[:,1])

# Signed heat
curves = [
					 [809, 808, 830, 783, 843, 181, 696, 781, 302, 937, 793, 26, 1025, 704, 1047, 1088, 487, 981, 1239, 1238, 402, 57, 1189, 161, 471, 529, 560, 316, 766, 832]
				 ]
signed_dist = solver.compute_signed_distance(curves, basisN, 
	preserve_source_normals=False, level_set_constraint="ZeroSet", soft_level_set_weight=-1)
ps_cloud.add_scalar_quantity("signed distance", signed_dist)

ps.show()
