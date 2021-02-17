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
