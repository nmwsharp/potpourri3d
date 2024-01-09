#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Eigen/Dense"

namespace py = pybind11;

using namespace geometrycentral;
using namespace geometrycentral::surface;


// For overloaded functions, with C++11 compiler only
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;


// A wrapper class for the heat method solver, which exposes Eigen in/out
class HeatMethodDistanceEigen {

public:
  HeatMethodDistanceEigen(DenseMatrix<double> verts, DenseMatrix<int64_t> faces, double tCoef = 1.0,
                          bool useRobustLaplacian = true) {

    // Construct the internal mesh and geometry
    mesh.reset(new SurfaceMesh(faces));
    geom.reset(new VertexPositionGeometry(*mesh));
    for (size_t i = 0; i < mesh->nVertices(); i++) {
      for (size_t j = 0; j < 3; j++) {
        geom->inputVertexPositions[i][j] = verts(i, j);
      }
    }

    // Build the solver
    solver.reset(new HeatMethodDistanceSolver(*geom, tCoef, useRobustLaplacian));
  }

  // Solve for distance from a single vertex
  Vector<double> compute_distance(int64_t sourceVert) {
    VertexData<double> dist = solver->computeDistance(mesh->vertex(sourceVert));
    return dist.toVector();
  }

  // Solve for distance from a collection of vertices
  Vector<double> compute_distance_multisource(Vector<int64_t> sourceVerts) {
    std::vector<Vertex> sources;
    for (size_t i = 0; i < sourceVerts.rows(); i++) {
      sources.push_back(mesh->vertex(sourceVerts(i)));
    }
    VertexData<double> dist = solver->computeDistance(sources);
    return dist.toVector();
  }

private:
  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::unique_ptr<HeatMethodDistanceSolver> solver;
};


// A wrapper class for the vector heat method solver, which exposes Eigen in/out
class VectorHeatMethodEigen {

  // TODO use intrinsic triangulations here

public:
  VectorHeatMethodEigen(DenseMatrix<double> verts, DenseMatrix<int64_t> faces, double tCoef = 1.0) {

    // Construct the internal mesh and geometry
    mesh.reset(new ManifoldSurfaceMesh(faces));
    geom.reset(new VertexPositionGeometry(*mesh));
    for (size_t i = 0; i < mesh->nVertices(); i++) {
      for (size_t j = 0; j < 3; j++) {
        geom->inputVertexPositions[i][j] = verts(i, j);
      }
    }

    // Build the solver
    solver.reset(new VectorHeatMethodSolver(*geom, tCoef));
  }

  // Extend scalars from a collection of vertices
  Vector<double> extend_scalar(Vector<int64_t> sourceVerts, Vector<double> values) {
    std::vector<std::tuple<Vertex, double>> sources;
    for (size_t i = 0; i < sourceVerts.rows(); i++) {
      sources.emplace_back(mesh->vertex(sourceVerts(i)), values(i));
    }
    VertexData<double> ext = solver->extendScalar(sources);
    return ext.toVector();
  }


  // Returns an extrinsic representation of the tangent frame being used internally, as X/Y/N vectors.
  std::tuple<DenseMatrix<double>, DenseMatrix<double>, DenseMatrix<double>> get_tangent_frames() {

    // Just in case we don't already have it
    geom->requireVertexNormals();
    geom->requireVertexTangentBasis();

    // unpack
    VertexData<Vector3> basisX(*mesh);
    VertexData<Vector3> basisY(*mesh);
    for (Vertex v : mesh->vertices()) {
      basisX[v] = geom->vertexTangentBasis[v][0];
      basisY[v] = geom->vertexTangentBasis[v][1];
    }

    return std::tuple<DenseMatrix<double>, DenseMatrix<double>, DenseMatrix<double>>(
        EigenMap<double, 3>(basisX), EigenMap<double, 3>(basisY), EigenMap<double, 3>(geom->vertexNormals));
  }

  SparseMatrix<std::complex<double>> get_connection_laplacian() {
    geom->requireVertexConnectionLaplacian();
    SparseMatrix<std::complex<double>> Lconn = geom->vertexConnectionLaplacian;
    geom->unrequireVertexConnectionLaplacian();
    return Lconn;
  }

  // TODO think about how to pass tangent frames around
  DenseMatrix<double> transport_tangent_vectors(Vector<int64_t> sourceVerts, DenseMatrix<double> values) {

    // Pack it as a Vector2
    std::vector<std::tuple<Vertex, Vector2>> sources;
    for (size_t i = 0; i < sourceVerts.rows(); i++) {
      sources.emplace_back(mesh->vertex(sourceVerts(i)), Vector2{values(i, 0), values(i, 1)});
    }
    VertexData<Vector2> ext = solver->transportTangentVectors(sources);

    return EigenMap<double, 2>(ext);
  }

  DenseMatrix<double> transport_tangent_vector(int64_t sourceVert, DenseMatrix<double> values) {

    // Pack it as a Vector2
    std::vector<std::tuple<Vertex, Vector2>> sources;
    sources.emplace_back(mesh->vertex(sourceVert), Vector2{values(0), values(1)});
    VertexData<Vector2> ext = solver->transportTangentVectors(sources);

    return EigenMap<double, 2>(ext);
  }


  DenseMatrix<double> compute_log_map(int64_t sourceVert) {
    return EigenMap<double, 2>(solver->computeLogMap(mesh->vertex(sourceVert)));
  }

private:
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::unique_ptr<VectorHeatMethodSolver> solver;
};

// A wrapper class for flip-based geodesics
class EdgeFlipGeodesicsManager {

public:
  EdgeFlipGeodesicsManager(DenseMatrix<double> verts, DenseMatrix<int64_t> faces) {

    // Construct the internal mesh and geometry
    mesh.reset(new ManifoldSurfaceMesh(faces));
    geom.reset(new VertexPositionGeometry(*mesh));
    for (size_t i = 0; i < mesh->nVertices(); i++) {
      for (size_t j = 0; j < 3; j++) {
        geom->inputVertexPositions[i][j] = verts(i, j);
      }
    }

    // Build the solver
    flipNetwork.reset(new FlipEdgeNetwork(*mesh, *geom, {}));
    flipNetwork->posGeom = geom.get();
    flipNetwork->supportRewinding = true;
  }

  // Generate a point-to-point geodesic by straightening a Dijkstra path
  DenseMatrix<double> find_geodesic_path(int64_t startVert, int64_t endVert) {

    // Get an initial dijkstra path
    std::vector<Halfedge> dijkstraPath = shortestEdgePath(*geom, mesh->vertex(startVert), mesh->vertex(endVert));

    if (startVert == endVert) {
      throw std::runtime_error("start and end vert are same");
    }
    if (dijkstraPath.empty()) {
      throw std::runtime_error("vertices lie on disconnected components of the surface");
    }

    // Reinitialize the ede network to contain this path
    flipNetwork->reinitializePath({dijkstraPath});

    // Straighten the path to geodesic
    flipNetwork->iterativeShorten();

    // Extract the path and store it in the vector
    std::vector<Vector3> path3D = flipNetwork->getPathPolyline3D().front();
    DenseMatrix<double> out(path3D.size(), 3);
    for (size_t i = 0; i < path3D.size(); i++) {
      for (size_t j = 0; j < 3; j++) {
        out(i, j) = path3D[i][j];
      }
    }

    // Be kind, rewind
    flipNetwork->rewind();

    return out;
  }

  // Generate a point-to-point geodesic by straightening a poly-geodesic path
  DenseMatrix<double> find_geodesic_path_poly(std::vector<int64_t> verts) {

    // Convert to a list of vertices
    std::vector<Halfedge> halfedges;

    for (size_t i = 0; i + 1 < verts.size(); i++) {
      Vertex vA = mesh->vertex(verts[i]);
      Vertex vB = mesh->vertex(verts[i + 1]);
      std::vector<Halfedge> dijkstraPath = shortestEdgePath(*geom, vA, vB);

      // validate
      if (vA == vB) {
        throw std::runtime_error("consecutive vertices are same");
      }
      if (dijkstraPath.empty()) {
        throw std::runtime_error("vertices lie on disconnected components of the surface");
      }

      halfedges.insert(halfedges.end(), dijkstraPath.begin(), dijkstraPath.end());
    }

    // Reinitialize the ede network to contain this path
    flipNetwork->reinitializePath({halfedges});

    // Straighten the path to geodesic
    flipNetwork->iterativeShorten();

    // Extract the path and store it in the vector
    std::vector<Vector3> path3D = flipNetwork->getPathPolyline3D().front();
    DenseMatrix<double> out(path3D.size(), 3);
    for (size_t i = 0; i < path3D.size(); i++) {
      for (size_t j = 0; j < 3; j++) {
        out(i, j) = path3D[i][j];
      }
    }

    // Be kind, rewind
    flipNetwork->rewind();

    return out;
  }


  // Generate a point-to-point geodesic loop by straightening a poly-geodesic path
  DenseMatrix<double> find_geodesic_loop(std::vector<int64_t> verts) {

    // Convert to a list of vertices
    std::vector<Halfedge> halfedges;

    for (size_t i = 0; i < verts.size(); i++) {
      Vertex vA = mesh->vertex(verts[i]);
      Vertex vB = mesh->vertex(verts[(i + 1) % verts.size()]);
      std::vector<Halfedge> dijkstraPath = shortestEdgePath(*geom, vA, vB);

      // validate
      if (vA == vB) {
        throw std::runtime_error("consecutive vertices are same");
      }
      if (dijkstraPath.empty()) {
        throw std::runtime_error("vertices lie on disconnected components of the surface");
      }

      halfedges.insert(halfedges.end(), dijkstraPath.begin(), dijkstraPath.end());
    }

    // Reinitialize the ede network to contain this path
    flipNetwork->reinitializePath({halfedges});

    // Straighten the path to geodesic
    flipNetwork->iterativeShorten();

    // Extract the path and store it in the vector
    std::vector<Vector3> path3D = flipNetwork->getPathPolyline3D().front();
    DenseMatrix<double> out(path3D.size(), 3);
    for (size_t i = 0; i < path3D.size(); i++) {
      for (size_t j = 0; j < 3; j++) {
        out(i, j) = path3D[i][j];
      }
    }

    // Be kind, rewind
    flipNetwork->rewind();

    return out;
  }

  // Generate a geodesic by tracing from a vertex along a tangent direction
  DenseMatrix<double> trace_geodesic(int64_t startVert, double direction_x, double direction_y,
                                     size_t max_iters = INVALID_IND) {
    TraceGeodesicResult result =
        traceGeodesic(*geom, SurfacePoint(mesh->vertex(startVert)), Vector2{direction_x, direction_y},
                      TraceOptions{true, false, NULL, max_iters});

    if (!result.hasPath) {
      throw std::runtime_error("geodesic trace encountered an error");
    }

    // Extract the path and store it in the vector
    DenseMatrix<double> out(result.pathPoints.size(), 3);
    for (size_t i = 0; i < result.pathPoints.size(); i++) {
      Vector3 point = result.pathPoints[i].interpolate(geom->vertexPositions);
      for (size_t j = 0; j < 3; j++) {
        out(i, j) = point[j];
      }
    }

    return out;
  }

private:
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::unique_ptr<FlipEdgeNetwork> flipNetwork;
};


// Actual binding code
// clang-format off
void bind_mesh(py::module& m) {

  py::class_<HeatMethodDistanceEigen>(m, "MeshHeatMethodDistance")
        .def(py::init<DenseMatrix<double>, DenseMatrix<int64_t>, double, bool>())
        .def("compute_distance", &HeatMethodDistanceEigen::compute_distance, py::arg("source_vert"))
        .def("compute_distance_multisource", &HeatMethodDistanceEigen::compute_distance_multisource, py::arg("source_verts"));
 

  py::class_<VectorHeatMethodEigen>(m, "MeshVectorHeatMethod")
        .def(py::init<DenseMatrix<double>, DenseMatrix<int64_t>, double>())
        .def("extend_scalar", &VectorHeatMethodEigen::extend_scalar, py::arg("source_verts"), py::arg("values"))
        .def("get_tangent_frames", &VectorHeatMethodEigen::get_tangent_frames)
        .def("get_connection_laplacian", &VectorHeatMethodEigen::get_connection_laplacian)
        .def("transport_tangent_vector", &VectorHeatMethodEigen::transport_tangent_vector, py::arg("source_vert"), py::arg("vector"))
        .def("transport_tangent_vectors", &VectorHeatMethodEigen::transport_tangent_vectors, py::arg("source_verts"), py::arg("vectors"))
        .def("compute_log_map", &VectorHeatMethodEigen::compute_log_map, py::arg("source_vert"));


  py::class_<EdgeFlipGeodesicsManager>(m, "EdgeFlipGeodesicsManager")
        .def(py::init<DenseMatrix<double>, DenseMatrix<int64_t>>())
        .def("find_geodesic_path", &EdgeFlipGeodesicsManager::find_geodesic_path, py::arg("source_vert"), py::arg("target_vert"))
        .def("find_geodesic_path_poly", &EdgeFlipGeodesicsManager::find_geodesic_path_poly, py::arg("vert_list"))
        .def("find_geodesic_loop", &EdgeFlipGeodesicsManager::find_geodesic_loop, py::arg("vert_list"))
        .def("trace_geodesic", &EdgeFlipGeodesicsManager::trace_geodesic, py::arg("start_vert"), py::arg("direction_x"), py::arg("direction_y"), py::arg("max_iters"));

  //m.def("read_mesh", &read_mesh, "Read a mesh from file.", py::arg("filename"));
}
