import numpy as np
import potpourri3d_bindings as pp3db

from .core import *

class PointCloudHeatSolver():

    def __init__(self, P, t_coef=1.):
        validate_points(P)
        self.bound_solver = pp3db.PointCloudHeatSolver(P, t_coef)

    def compute_distance(self, p_ind):
        return self.bound_solver.compute_distance(p_ind)
    
    def compute_distance_multisource(self, p_inds):
        return self.bound_solver.compute_distance_multisource(p_inds)
    
    def extend_scalar(self, p_inds, values):
        if len(p_inds) != len(values):
            raise ValueError("source point indices and values array should be same shape")
        return self.bound_solver.extend_scalar(p_inds, values)
    
    def get_tangent_frames(self):
        return self.bound_solver.get_tangent_frames()
    
    def transport_tangent_vector(self, p_ind, vector):
        if len(vector) != 2:
            raise ValueError("vector should be a 2D tangent vector")
        return self.bound_solver.transport_tangent_vector(p_ind, vector)
    
    def transport_tangent_vectors(self, p_inds, vectors):
        if len(p_inds) != len(vectors):
            raise ValueError("source point indices and values array should be same length")
        return self.bound_solver.transport_tangent_vectors(p_inds, vectors)
    
    def compute_log_map(self, p_ind):
        return self.bound_solver.compute_log_map(p_ind)
