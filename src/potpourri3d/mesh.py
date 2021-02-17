import numpy as np
import potpourri3d_bindings as pp3db

from .core import *

class MeshHeatMethodDistanceSolver():

    def __init__(self, V, F, t_coef=1., use_robust=True):
        validate_mesh(V, F, force_triangular=True)
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

    def __init__(self, V, F, t_coef=1.):
        validate_mesh(V, F, force_triangular=True)
        self.bound_solver = pp3db.MeshVectorHeatMethod(V, F, t_coef)

    def extend_scalar(self, v_inds, values):
        if len(v_inds) != len(values):
            raise ValueError("source vertex indices and values array should be same shape")
        return self.bound_solver.extend_scalar(v_inds, values)
    
    def get_tangent_frames(self):
        return self.bound_solver.get_tangent_frames()
    
    def transport_tangent_vector(self, v_ind, vector):
        if len(vector) != 2:
            raise ValueError("vector should be a 2D tangent vector")
        return self.bound_solver.transport_tangent_vector(v_ind, vector)
    
    def transport_tangent_vectors(self, v_inds, vectors):
        if len(v_inds) != len(vectors):
            raise ValueError("source vertex indices and values array should be same length")
        return self.bound_solver.transport_tangent_vectors(v_inds, vectors)
    
    def compute_log_map(self, v_ind):
        return self.bound_solver.compute_log_map(v_ind)
