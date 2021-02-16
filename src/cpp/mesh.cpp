#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

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
    std::vector<std::tuple<Vertex,double>> sources;
    for (size_t i = 0; i < sourceVerts.rows(); i++) {
      sources.emplace_back(mesh->vertex(sourceVerts(i)), values(i));
    }
    VertexData<double> ext = solver->extendScalar(sources);
    return ext.toVector();
  }

  /* TODO think about how to pass tangent frames around
  DenseMatrix<double> transport_tangent_vectors(Vector<int64_t> sourceVerts, DenseMatrix<double> values);
  DenseMatrix<double> compute_log_map(int64_t sourceVert);
  */

private:
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::unique_ptr<VectorHeatMethodSolver> solver;
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
        .def("extend_scalar", &VectorHeatMethodEigen::extend_scalar, py::arg("source_verts"), py::arg("values"));


  //m.def("read_mesh", &read_mesh, "Read a mesh from file.", py::arg("filename"));
}