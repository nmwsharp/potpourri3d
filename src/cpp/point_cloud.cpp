#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/pointcloud/point_cloud_heat_solver.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "Eigen/Dense"

namespace py = pybind11;

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace geometrycentral::pointcloud;


// For overloaded functions, with C++11 compiler only
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;


// A wrapper class for the heat method solver, which exposes Eigen in/out
class PointCloudHeatSolverEigen {

public:
  PointCloudHeatSolverEigen(DenseMatrix<double> points, double tCoef = 1.0) {

    // Construct the internal cloud and geometry
    cloud.reset(new PointCloud(points.rows()));
    geom.reset(new PointPositionGeometry(*cloud));
    for (size_t i = 0; i < cloud->nPoints(); i++) {
      for (size_t j = 0; j < 3; j++) {
        geom->positions[i][j] = points(i, j);
      }
    }

    // Build the solver
    solver.reset(new PointCloudHeatSolver(*cloud, *geom, tCoef));
  }

  // Solve for distance from a single point
  Vector<double> compute_distance(int64_t sourcePoint) {
    PointData<double> dist = solver->computeDistance(cloud->point(sourcePoint));
    return dist.toVector();
  }

  // Solve for distance from a collection of points
  Vector<double> compute_distance_multisource(Vector<int64_t> sourcePoints) {
    std::vector<Point> sources;
    for (size_t i = 0; i < sourcePoints.rows(); i++) {
      sources.push_back(cloud->point(sourcePoints(i)));
    }
    PointData<double> dist = solver->computeDistance(sources);
    return dist.toVector();
  }


  Vector<double> extend_scalar(Vector<int64_t> sourcePoints, Vector<double> values) {
    std::vector<std::tuple<Point, double>> sources;
    for (size_t i = 0; i < sourcePoints.rows(); i++) {
      sources.emplace_back(cloud->point(sourcePoints(i)), values(i));
    }
    PointData<double> ext = solver->extendScalars(sources);
    return ext.toVector();
  }
  
  // Returns an extrinsic representation of the tangent frame being used internally, as X/Y/N vectors.
  std::tuple<DenseMatrix<double>, DenseMatrix<double>, DenseMatrix<double>> get_tangent_frames() {

    // Just in case we don't already have it
    geom->requireNormals();
    geom->requireTangentBasis();

    // unpack
    PointData<Vector3> basisX(*cloud);
    PointData<Vector3> basisY(*cloud);
    for (Point v : cloud->points()) {
      basisX[v] = geom->tangentBasis[v][0];
      basisY[v] = geom->tangentBasis[v][1];
    }

    return std::tuple<DenseMatrix<double>, DenseMatrix<double>, DenseMatrix<double>>(
        EigenMap<double, 3>(basisX), EigenMap<double, 3>(basisY), EigenMap<double, 3>(geom->normals));
  }

  DenseMatrix<double> transport_tangent_vectors(Vector<int64_t> sourcePoints, DenseMatrix<double> values) {

    // Pack it as a Vector2
    std::vector<std::tuple<Point, Vector2>> sources;
    for (size_t i = 0; i < sourcePoints.rows(); i++) {
      sources.emplace_back(cloud->point(sourcePoints(i)), Vector2{values(i, 0), values(i, 1)});
    }
    PointData<Vector2> ext = solver->transportTangentVectors(sources);

    return EigenMap<double, 2>(ext);
  }

  DenseMatrix<double> transport_tangent_vector(int64_t sourcePoint, DenseMatrix<double> values) {

    // Pack it as a Vector2
    std::vector<std::tuple<Point, Vector2>> sources;
    sources.emplace_back(cloud->point(sourcePoint), Vector2{values(0), values(1)});
    PointData<Vector2> ext = solver->transportTangentVectors(sources);

    return EigenMap<double, 2>(ext);
  }


  DenseMatrix<double> compute_log_map(int64_t sourcePoint) {
    return EigenMap<double, 2>(solver->computeLogMap(cloud->point(sourcePoint)));
  }

private:
  std::unique_ptr<PointCloud> cloud;
  std::unique_ptr<PointPositionGeometry> geom;
  std::unique_ptr<PointCloudHeatSolver> solver;
};


// Actual binding code
// clang-format off
void bind_point_cloud(py::module& m) {
  
  py::class_<PointCloudHeatSolverEigen>(m, "PointCloudHeatSolver")
        .def(py::init<DenseMatrix<double>, double>())
        .def("compute_distance", &PointCloudHeatSolverEigen::compute_distance, py::arg("source_point"))
        .def("compute_distance_multisource", &PointCloudHeatSolverEigen::compute_distance_multisource, py::arg("source_points"))
        .def("extend_scalar", &PointCloudHeatSolverEigen::extend_scalar, py::arg("source_points"), py::arg("source_values"))
        .def("get_tangent_frames", &PointCloudHeatSolverEigen::get_tangent_frames)
        .def("transport_tangent_vector", &PointCloudHeatSolverEigen::transport_tangent_vector, py::arg("source_point"), py::arg("vector"))
        .def("transport_tangent_vectors", &PointCloudHeatSolverEigen::transport_tangent_vectors, py::arg("source_points"), py::arg("vectors"))
        .def("compute_log_map", &PointCloudHeatSolverEigen::compute_log_map, py::arg("source_point"));
}
