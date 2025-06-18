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

    def compute_signed_distance(self, curves, cloud_normals, 
                                preserve_source_normals=False, level_set_constraint="ZeroSet", soft_level_set_weight=-1):
        return self.bound_solver.compute_signed_distance(curves, cloud_normals, 
            preserve_source_normals, level_set_constraint, soft_level_set_weight)


class PointCloudLocalTriangulation():

    def __init__(self, P, with_degeneracy_heuristic=True):
        validate_points(P)
        self.bound_triangulation = pp3db.PointCloudLocalTriangulation(P, with_degeneracy_heuristic)

    def get_local_triangulation(self):
        """Return the local point cloud triangulation
        
        The out matrix has the following convention:
            size: num_points, max_neighs, 3. max_neighs is the maximum number of neighbors
            out[point_idx, neigh_idx, :] are the indices of the 3 neighbors
            -1 is used as the fill value for unused elements if num_neighs < max_neighs for a point
        """
        out = self.bound_triangulation.get_local_triangulation()
        assert out.shape[-1] % 3 == 0
        max_neighs = out.shape[-1] // 3
        out = np.reshape(out, [-1, max_neighs, 3])
        return out
