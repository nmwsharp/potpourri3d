import numpy as np
import potpourri3d_bindings as pp3db

import scipy
import scipy.sparse

from .core import *

class MeshFastMarchingDistanceSolver():

    def __init__(self, V, F):
        validate_mesh(V, F, force_triangular=True, test_indices=True)
        self.bound_solver = pp3db.MeshFastMarchingDistance(V, F)

    def compute_distance(self, curves, distances=[], sign=False): 
        return self.bound_solver.compute_distance(curves, distances, sign)

class MeshMarchingTrianglesSolver():

    def __init__(self, V, F):
        validate_mesh(V, F, force_triangular=True, test_indices=True)
        self.bound_solver = pp3db.MeshMarchingTriangles(V, F)

    def marching_triangles(self, u, isoval=0.): 
        return self.bound_solver.marching_triangles(u, isoval)

def marching_triangles(V, F, u, isoval=0.): 
    solver = MeshMarchingTrianglesSolver(V, F)
    return solver.marching_triangles(u, isoval)

class MeshHeatMethodDistanceSolver():

    def __init__(self, V, F, t_coef=1., use_robust=True):
        validate_mesh(V, F, force_triangular=True, test_indices=True)
        self.bound_solver = pp3db.MeshHeatMethodDistance(V, F, t_coef, use_robust)

    def compute_distance(self, v_ind):
        return self.bound_solver.compute_distance(v_ind)
    
    def compute_distance_multisource(self, v_inds):
        return self.bound_solver.compute_distance_multisource(v_inds)

def compute_distance(V, F, v_ind): 
    solver = MeshHeatMethodDistanceSolver(V, F)
    return solver.compute_distance(v_ind)

def compute_distance_multisource(V, F, v_inds):
    solver = MeshHeatMethodDistanceSolver(V, F)
    return solver.compute_distance_multisource(v_inds)


class MeshVectorHeatSolver():

    def __init__(self, V, F, t_coef=1., use_intrinsic_delaunay=True):
        validate_mesh(V, F, force_triangular=True, test_indices=True)
        self.bound_solver = pp3db.MeshVectorHeatMethod(V, F, t_coef, use_intrinsic_delaunay)

    def extend_scalar(self, v_inds, values):
        if len(v_inds) != len(values):
            raise ValueError("source vertex indices and values array should be same shape")
        return self.bound_solver.extend_scalar(v_inds, values)
    
    def get_tangent_frames(self):
        return self.bound_solver.get_tangent_frames()
    
    def get_connection_laplacian(self):
        return self.bound_solver.get_connection_laplacian()
    
    def transport_tangent_vector(self, v_ind, vector):
        if len(vector) != 2:
            raise ValueError("vector should be a 2D tangent vector")
        return self.bound_solver.transport_tangent_vector(v_ind, vector)
    
    def transport_tangent_vectors(self, v_inds, vectors):
        if len(v_inds) != len(vectors):
            raise ValueError("source vertex indices and values array should be same length")
        return self.bound_solver.transport_tangent_vectors(v_inds, vectors)
    
    def compute_log_map(self, v_ind, strategy="AffineLocal"):
        return self.bound_solver.compute_log_map(v_ind, strategy)


class MeshSignedHeatSolver():

    def __init__(self, V, F, t_coef=1.):
        validate_mesh(V, F, force_triangular=True, test_indices=True)
        self.bound_solver = pp3db.MeshSignedHeatMethod(V, F, t_coef)

    def compute_distance(self, curves, curve_signs=[], points=[],
                         preserve_source_normals=False, level_set_constraint="ZeroSet", soft_level_set_weight=-1):
        
        '''
        Args:
            curves [list]: List of curves. Each curve is a list of tuples (element_indices, barycentric_coords).
            curve_signs [list]: list of bools, indicating whether each curve in `curves' 
                                is signed (true) or unsigned (false). 
                                By default, curves are assumed to be signed (oriented).
        '''
        return self.bound_solver.compute_distance(curves, curve_signs, points, 
            preserve_source_normals, level_set_constraint, soft_level_set_weight)


class PolygonMeshHeatSolver():

    def __init__(self, V, F, t_coef=1.):
        validate_mesh(V, F, force_triangular=False, test_indices=True)
        self.bound_solver = pp3db.PolygonMeshHeatSolver(V, F, t_coef)

    def compute_distance(self, v_inds):
        return self.bound_solver.compute_distance(v_inds)

    def extend_scalar(self, v_inds, values):
        if len(v_inds) != len(values):
            raise ValueError("source vertex indices and values array should be same shape")
        return self.bound_solver.extend_scalar(v_inds, values)

    def get_tangent_frames(self):
        return self.bound_solver.get_tangent_frames()
    
    def transport_tangent_vectors(self, v_inds, vectors):
        if len(v_inds) != len(vectors):
            raise ValueError("source vertex indices and values array should be same length")
        return self.bound_solver.transport_tangent_vectors(v_inds, vectors)

    def compute_signed_distance(self, curves, level_set_constraint="ZeroSet"):
        return self.bound_solver.compute_signed_distance(curves, level_set_constraint)


class EdgeFlipGeodesicSolver():

    def __init__(self, V, F, t_coef=1.):
        validate_mesh(V, F, force_triangular=True)
        self.bound_solver = pp3db.EdgeFlipGeodesicsManager(V, F)

    def find_geodesic_path(self, v_start, v_end, max_iterations=None, max_relative_length_decrease=None):
        if max_iterations is None:
            max_iterations = 2**63-1
        if max_relative_length_decrease is None:
            max_relative_length_decrease = 0.

        return self.bound_solver.find_geodesic_path(v_start, v_end, max_iterations, max_relative_length_decrease)
    
    def find_geodesic_path_poly(self, v_list, max_iterations=None, max_relative_length_decrease=None):
        if max_iterations is None:
            max_iterations = 2**63-1
        if max_relative_length_decrease is None:
            max_relative_length_decrease = 0.

        return self.bound_solver.find_geodesic_path_poly(v_list, max_iterations, max_relative_length_decrease)
    
    def find_geodesic_loop(self, v_list, max_iterations=None, max_relative_length_decrease=None):
        if max_iterations is None:
            max_iterations = 2**63-1
        if max_relative_length_decrease is None:
            max_relative_length_decrease = 0.

        return self.bound_solver.find_geodesic_loop(v_list, max_iterations, max_relative_length_decrease)

class GeodesicTracer():

    def __init__(self, V, F, t_coef=1.):
        validate_mesh(V, F, force_triangular=True)
        self.bound_tracer = pp3db.GeodesicTracer(V, F)

    def trace_geodesic_from_vertex(self, start_vert, direction_xyz, max_iterations=None):
        if max_iterations is None:
            max_iterations = 2**63-1

        direction_xyz = np.array(direction_xyz)

        return self.bound_tracer.trace_geodesic_from_vertex(start_vert, direction_xyz, max_iterations)

    def trace_geodesic_from_face(self, start_face, bary_coords, direction_xyz, max_iterations=None):
        if max_iterations is None:
            max_iterations = 2**63-1

        bary_coords = np.array(bary_coords)
        direction_xyz = np.array(direction_xyz)

        return self.bound_tracer.trace_geodesic_from_face(start_face, bary_coords, direction_xyz, max_iterations)
    

def cotan_laplacian(V, F, denom_eps=0.):
    validate_mesh(V, F, force_triangular=True)
    nV = V.shape[0]

    mat_i = []
    mat_j = []
    mat_data = []
    for i in range(3):
      
        # Gather indices and compute cotan weight (via dot() / cross() formula)
        inds_i = F[:,i]
        inds_j = F[:,(i+1)%3]
        inds_k = F[:,(i+2)%3]
        vec_ki = V[inds_i,:] - V[inds_k,:]
        vec_kj = V[inds_j,:] - V[inds_k,:]
        dots = np.sum(vec_ki * vec_kj, axis=1)
        cross_mags =  np.linalg.norm(np.cross(vec_ki, vec_kj), axis=1)
        cotans = 0.5 * dots / (cross_mags + denom_eps)

        # Add the four matrix entries from this weight

        mat_i.append(inds_i)
        mat_j.append(inds_i)
        mat_data.append(cotans)
        
        mat_i.append(inds_j)
        mat_j.append(inds_j)
        mat_data.append(cotans)
        
        mat_i.append(inds_i)
        mat_j.append(inds_j)
        mat_data.append(-cotans)
        
        mat_i.append(inds_j)
        mat_j.append(inds_i)
        mat_data.append(-cotans)

    # Concatenate the arrays to single lists
    mat_i = np.concatenate(mat_i)
    mat_j = np.concatenate(mat_j)
    mat_data = np.concatenate(mat_data)

    L_coo = scipy.sparse.coo_matrix((mat_data, (mat_i, mat_j)), shape=(nV, nV))

    return L_coo.tocsr()

def face_areas(V, F):
    validate_mesh(V, F, force_triangular=True)

    vec_ij = V[F[:,1],:] - V[F[:,0],:]
    vec_ik = V[F[:,2],:] - V[F[:,0],:]

    areas = 0.5 * np.linalg.norm(np.cross(vec_ij, vec_ik), axis=1)
    
    return areas

def vertex_areas(V, F):
    validate_mesh(V, F, force_triangular=True)
    nV = V.shape[0]

    face_area = face_areas(V, F)

    vertex_area = np.zeros(V.shape[0])
    for i in range(3):
        vertex_area += np.bincount(F[:,i], face_area, minlength=nV)
    vertex_area /= 3.
    
    return vertex_area

def edges(V, F):
    validate_mesh(V, F, force_triangular=False)
    return pp3db.edges(V, F)
