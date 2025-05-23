import numpy as np
import potpourri3d_bindings as pp3db

from .core import *

def read_mesh(filename):
    V, F = pp3db.read_mesh(filename)
    V = np.ascontiguousarray(V)
    F = np.ascontiguousarray(F)
    return V, F

def read_polygon_mesh(filename):
    V, F = pp3db.read_polygon_mesh(filename)
    V = np.ascontiguousarray(V)
    return V, F

def write_mesh(V, F, filename, UV_coords=None, UV_type=None):
    # TODO generalize this to take indexed UVs
    # (the underlying geometry-central writer needs to support it first)

    validate_mesh(V, F, test_indices=True)

    if UV_type is None:

        pp3db.write_mesh(V, F, filename)

    elif UV_type == 'per-vertex':

        if len(UV_coords.shape) != 2 or UV_coords.shape[0] != V.shape[0] or UV_coords.shape[1] != 2:
            raise ValueError("UV_coords should be a 2d Vx2 numpy array")

        pp3db.write_mesh_pervertex_uv(V, F, UV_coords, filename)

    elif UV_type == 'per-face':

        if len(UV_coords.shape) != 2 or UV_coords.shape[0] != F.shape[0] or UV_coords.shape[1] != 2:
            raise ValueError("UV_coords should be a 2d Fx2 numpy array")

        pp3db.write_mesh_perface_uv(V, F, UV_coords, filename)

    elif UV_type == 'per-corner':

        if len(UV_coords.shape) != 2 or UV_coords.shape[0] != F.shape[0]*F.shape[1] or UV_coords.shape[1] != 2:
            raise ValueError("UV_coords should be a 2d Fx2 numpy array")

        pp3db.write_mesh_percorner_uv(V, F, UV_coords, filename)

    else:
        raise ValueError(f"unrecognized value for UV_type: {UV_type}. Should be one of: [None, 'per-vertex', 'per-face', 'per-corner']")


def read_point_cloud(filename):
    V = pp3db.read_point_cloud(filename)
    V = np.ascontiguousarray(V)
    return V

def write_point_cloud(V, filename):
    validate_points(V)
    pp3db.write_point_cloud(V, filename)
