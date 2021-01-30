import os, sys

import polyscope as ps
import numpy as np
# import scipy.sparse.linalg as sla

# Path to where the bindings live
sys.path.append(os.path.join(os.path.dirname(__file__), "../build/"))
sys.path.append(os.path.join(os.path.dirname(__file__), "../src/"))

import potpourri3d as pp3d

ps.init()

# Read input
V, F = pp3d.read_mesh("/Users/nick/mesh/spot.obj")

ps_mesh = ps.register_surface_mesh("mesh", V, F)

ps.show()
