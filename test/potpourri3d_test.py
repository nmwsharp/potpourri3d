import unittest
import os
import sys
import os.path as path 
import numpy as np
import scipy

# Path to where the bindings live
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
if os.name == 'nt': # if Windows
    # handle default location where VS puts binary
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build", "Debug")))
else:
    # normal / unix case
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "build")))


import potpourri3d as pp3d

asset_path = os.path.abspath(os.path.dirname(__file__))

def generate_verts(n_pts=999):
    np.random.seed(777)        
    return np.random.rand(n_pts, 3)

def generate_faces(n_pts=999):
    # n_pts should be a multiple of 3 for indexing to work out
    np.random.seed(777)        
    rand_faces = np.random.randint(0, n_pts, size=(2*n_pts,3))
    coverage_faces = np.arange(n_pts).reshape(-1, 3)
    faces = np.vstack((rand_faces, coverage_faces))
    return faces

def is_symmetric(A, eps=1e-6):
    resid = A - A.T
    return np.all(np.abs(resid.data) < eps)

def is_nonnegative(A, eps=1e-6):
    return np.all(A.data > -eps)

class TestCore(unittest.TestCase):
    
    def test_write_read_mesh(self):

        for ext in ['obj']:

            V = generate_verts()
            F = generate_faces()

            fname = "test." + ext

            # write
            pp3d.write_mesh(V,F,fname)

            Vnew, Fnew = pp3d.read_mesh(fname)
            
            self.assertLess(np.amax(np.abs(V-Vnew)), 1e-6)
            self.assertTrue((F==Fnew).all())


            # smoke test various UV writers
            UV_vert = V[:,:2]
            pp3d.write_mesh(V,F,fname,UV_coords=UV_vert, UV_type='per-vertex')

            UV_face = F[:,:2] * .3
            pp3d.write_mesh(V,F,fname,UV_coords=UV_face, UV_type='per-face')
            
            UV_corner = np.zeros((F.shape[0]*F.shape[1],2))
            pp3d.write_mesh(V,F,fname,UV_coords=UV_corner, UV_type='per-corner')
    
    def test_write_read_point_cloud(self):

        for ext in ['obj', 'ply']:

            V = generate_verts()

            fname = "test_cloud." + ext

            # write
            pp3d.write_point_cloud(V, fname)

            Vnew = pp3d.read_point_cloud(fname)
            
            self.assertLess(np.amax(np.abs(V-Vnew)), 1e-6)

        # self.assertTrue(is_nonnegative(off_L)) # positive edge weights
        # self.assertGreater(L.sum(), -1e-5)
        # self.assertEqual(M.sum(), M.diagonal().sum())
    

    def test_mesh_heat_distance(self):

        V = generate_verts()
        F = generate_faces()
         
        # Test stateful version
        
        solver = pp3d.MeshHeatMethodDistanceSolver(V,F)
        dist = solver.compute_distance(7)
        self.assertEqual(dist.shape[0], V.shape[0])

        dist = solver.compute_distance_multisource([1,2,3])
        self.assertEqual(dist.shape[0], V.shape[0])
        
        
        # = Test one-off versions

        dist = pp3d.compute_distance(V,F,7)
        self.assertEqual(dist.shape[0], V.shape[0])

        dist = pp3d.compute_distance_multisource(V,F,[1,3,4])
        self.assertEqual(dist.shape[0], V.shape[0])
   

    def test_mesh_vector_heat(self):

        V, F = pp3d.read_mesh(os.path.join(asset_path, "bunny_small.ply"))

        for use_intrinsic_delaunay in [False, True]:

            solver = pp3d.MeshVectorHeatSolver(V, F, use_intrinsic_delaunay=use_intrinsic_delaunay)

            # Scalar extension
            ext = solver.extend_scalar([1, 22], [0., 6.])
            self.assertEqual(ext.shape[0], V.shape[0])
            self.assertGreaterEqual(np.amin(ext), 0.)
            
            # Get frames
            basisX, basisY, basisN = solver.get_tangent_frames()
            self.assertEqual(basisX.shape[0], V.shape[0])
            self.assertEqual(basisY.shape[0], V.shape[0])
            self.assertEqual(basisN.shape[0], V.shape[0])
            # TODO could check orthogonal
            
            # Get connection Laplacian
            L_conn = solver.get_connection_laplacian()
            self.assertTrue(isinstance(L_conn, scipy.sparse.csc_matrix))
            max_diag_imag = np.max(np.abs(L_conn.diagonal().imag))
            self.assertLess(max_diag_imag, 1e-4)

            # Vector heat (transport vector)
            ext = solver.transport_tangent_vector(1, [6., 6.])
            self.assertEqual(ext.shape[0], V.shape[0])
            self.assertEqual(ext.shape[1], 2)
            ext = solver.transport_tangent_vectors([1, 22], [[6., 6.], [3., 4.]])
            self.assertEqual(ext.shape[0], V.shape[0])
            self.assertEqual(ext.shape[1], 2)

            # Vector heat (log map)
            for strategy in [None, "VectorHeat", "AffineLocal", "AffineAdaptive"]:
                kwargs = {} if strategy is None else {"strategy": strategy}
                logmap = solver.compute_log_map(1, **kwargs)
                self.assertEqual(logmap.shape[0], V.shape[0])
                self.assertEqual(logmap.shape[1], 2)
    
    def test_mesh_cotan_laplace(self):

        V, F = pp3d.read_mesh(os.path.join(asset_path, "bunny_small.ply"))

        L = pp3d.cotan_laplacian(V,F) 

        self.assertEqual(L.shape[0],V.shape[0])
        self.assertEqual(L.shape[1],V.shape[0])

        self.assertLess(np.abs(np.sum(L)), 1e-6)
    
    def test_mesh_areas(self):

        V, F = pp3d.read_mesh(os.path.join(asset_path, "bunny_small.ply"))

        face_area = pp3d.face_areas(V,F)
        self.assertEqual(face_area.shape[0],F.shape[0])
        self.assertTrue(np.all(face_area >= 0))

        vert_area = pp3d.vertex_areas(V,F) 
        self.assertLess(np.abs(np.sum(face_area) - np.sum(vert_area)), 1e-6)

    def test_mesh_flip_geodesic(self):

        V, F = pp3d.read_mesh(os.path.join(asset_path, "bunny_small.ply"))
         
        # Test stateful version
        path_solver = pp3d.EdgeFlipGeodesicSolver(V,F)

        # Do a first path
        path_pts = path_solver.find_geodesic_path(v_start=14, v_end=22)
        self.assertEqual(len(path_pts.shape), 2)
        self.assertEqual(path_pts.shape[1], 3)
        path_pts = path_solver.find_geodesic_path(v_start=14, v_end=22, max_iterations=100, max_relative_length_decrease=0.5)

        # Do some more
        for i in range(5):
            path_pts = path_solver.find_geodesic_path(v_start=14, v_end=22+i)
            self.assertEqual(len(path_pts.shape), 2)
            self.assertEqual(path_pts.shape[1], 3)

        # Initialize with a compound path
        path_pts = path_solver.find_geodesic_path_poly([1173, 148, 870, 898])
        self.assertEqual(len(path_pts.shape), 2)
        self.assertEqual(path_pts.shape[1], 3)
        path_pts = path_solver.find_geodesic_path_poly([1173, 148, 870, 898], max_iterations=100, max_relative_length_decrease=0.5)

        # Do a loop
        loop_pts = path_solver.find_geodesic_loop([1173, 148, 870, 898])
        self.assertEqual(len(loop_pts.shape), 2)
        self.assertEqual(loop_pts.shape[1], 3)

        # Do another loop
        # this one contracts to a point
        loop_pts = path_solver.find_geodesic_loop([307, 757, 190]) 
        self.assertEqual(len(loop_pts.shape), 2)
        self.assertEqual(loop_pts.shape[1], 3)
        loop_pts = path_solver.find_geodesic_loop([307, 757, 190], max_iterations=100, max_relative_length_decrease=0.5)
    

    def test_geodesic_trace(self):

        V, F = pp3d.read_mesh(os.path.join(asset_path, "bunny_small.ply"))

        # Test stateful version
        tracer = pp3d.GeodesicTracer(V,F)

        # Trace from a vertex
        trace_pts = tracer.trace_geodesic_from_vertex(22, np.array((0.3, 0.5, 0.4)))
        self.assertEqual(len(trace_pts.shape), 2)
        self.assertEqual(trace_pts.shape[1], 3)

        trace_pts = tracer.trace_geodesic_from_vertex(22, np.array((0.3, 0.5, 0.4)), max_iterations=10)
        self.assertEqual(len(trace_pts.shape), 2)
        self.assertEqual(trace_pts.shape[1], 3)

        # Trace from a face
        trace_pts = tracer.trace_geodesic_from_face(31, np.array((0.1, 0.4, 0.5)), np.array((0.3, 0.5, 0.4)))
        self.assertEqual(len(trace_pts.shape), 2)
        self.assertEqual(trace_pts.shape[1], 3)

        trace_pts = tracer.trace_geodesic_from_face(31, np.array((0.1, 0.4, 0.5)), np.array((0.3, 0.5, 0.4)), max_iterations=10)
        self.assertEqual(len(trace_pts.shape), 2)
        self.assertEqual(trace_pts.shape[1], 3)

    def test_point_cloud_distance(self):

        P = generate_verts()
       
        solver = pp3d.PointCloudHeatSolver(P)
        
        dist = solver.compute_distance(7)
        self.assertEqual(dist.shape[0], P.shape[0])

        dist = solver.compute_distance_multisource([1,2,3])
        self.assertEqual(dist.shape[0], P.shape[0])

    def test_point_cloud_vector_heat(self):

        P = generate_verts()
       
        solver = pp3d.PointCloudHeatSolver(P)

        # Scalar extension
        ext = solver.extend_scalar([1, 22], [0., 6.])
        self.assertEqual(ext.shape[0], P.shape[0])
        self.assertGreaterEqual(np.amin(ext), 0.)
        
        # Get frames
        basisX, basisY, basisN = solver.get_tangent_frames()
        self.assertEqual(basisX.shape[0], P.shape[0])
        self.assertEqual(basisY.shape[0], P.shape[0])
        self.assertEqual(basisN.shape[0], P.shape[0])
        # TODO could check orthogonal

        # Vector heat (transport vector)
        ext = solver.transport_tangent_vector(1, [6., 6.])
        self.assertEqual(ext.shape[0], P.shape[0])
        self.assertEqual(ext.shape[1], 2)
        ext = solver.transport_tangent_vectors([1, 22], [[6., 6.], [3., 4.]])
        self.assertEqual(ext.shape[0], P.shape[0])
        self.assertEqual(ext.shape[1], 2)

        # Vector heat (log map)
        logmap = solver.compute_log_map(1)
        self.assertEqual(logmap.shape[0], P.shape[0])
        self.assertEqual(logmap.shape[1], 2)

    def test_point_cloud_local_triangulation(self):
        # Test local triangulation for a "cartwheel" pointcloud

        num = 31
        t = np.linspace(0, 2*np.pi, num-1, endpoint=False)
        points = np.concatenate([np.zeros([1, 3]), np.stack([np.cos(t), np.sin(t), 0*t], 1)], 0)

        pcl_local_tri = pp3d.PointCloudLocalTriangulation(points)
        idxs = pcl_local_tri.get_local_triangulation()

        def next_id(i):
            assert i != 0
            if i == num-1:
                return 1
            return i+1

        def prev_id(i):
            assert i != 0
            if i == 1:
                return num-1
            return i-1


        # Explicitly check for cartwheel
        # Use sets for order invariance

        res0 = set(tuple(r) for r in idxs[0])
        ref0 = set((0, j, next_id(j)) for j in range(1, num))
        self.assertEqual(res0, ref0)

        for i in range(1, num):
            self.assertTrue(np.all(idxs[i, 2:] == -1))
            
            res = set(tuple(r) for r in idxs[i, :2])
            ref = {(i, next_id(i), 0), (i, 0, prev_id(i))}
            self.assertEqual(res, ref)


if __name__ == '__main__':
    unittest.main()
