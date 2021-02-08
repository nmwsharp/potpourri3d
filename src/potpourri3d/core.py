import numpy as np

import potpourri3d_bindings as pp3db

def read_mesh(filename):

    ## Call the bindings
    V, F = pp3db.read_mesh(filename)

    ## Return the result
    return V, F

def write_mesh(V, F, filename):

    ## Call the bindings
    pp3db.write_mesh(V, F, filename)
