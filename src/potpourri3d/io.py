import numpy as np
import potpourri3d_bindings as pp3db

def read_mesh(filename):
    V, F = pp3db.read_mesh(filename)
    return V, F

def write_mesh(V, F, filename):
    pp3db.write_mesh(V, F, filename)

def read_point_cloud(filename):
    V = pp3db.read_point_cloud(filename)
    return V

def write_point_cloud(V, filename):
    pp3db.write_point_cloud(V, filename)
