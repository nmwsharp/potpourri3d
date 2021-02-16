#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
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

std::tuple<DenseMatrix<double>, DenseMatrix<int64_t>> read_mesh(std::string filename) {

  // Call the mesh reader
  SimplePolygonMesh pmesh(filename);

  if (pmesh.nFaces() == 0) throw std::runtime_error("read mesh has no faces");

  // Manually copy the vertex array
  DenseMatrix<double> V(pmesh.nVertices(), 3);
  for (size_t i = 0; i < pmesh.nVertices(); i++) {
    for (size_t j = 0; j < 3; j++) {
      V(i, j) = pmesh.vertexCoordinates[i][j];
    }
  }

  // Manually copy the face array
  size_t fDegree = pmesh.polygons[0].size();
  DenseMatrix<int64_t> F(pmesh.nFaces(), fDegree);
  for (size_t i = 0; i < pmesh.nFaces(); i++) {
    if (pmesh.polygons[i].size() != fDegree) throw std::runtime_error("read mesh faces are not all the same degree");
    for (size_t j = 0; j < fDegree; j++) {
      F(i, j) = pmesh.polygons[i][j];
    }
  }


  return std::make_tuple(V, F);
}

void write_mesh(DenseMatrix<double> verts, DenseMatrix<int64_t> faces, std::string filename) {

  // Copy in to the mesh object
  std::vector<Vector3> coords(verts.rows());
  for (size_t i = 0; i < verts.rows(); i++) {
    for (size_t j = 0; j < 3; j++) {
      coords[i][j] = verts(i, j);
    }
  }
  std::vector<std::vector<size_t>> polys(faces.rows());
  for (size_t i = 0; i < faces.rows(); i++) {
    polys[i].resize(faces.cols());
    for (size_t j = 0; j < faces.cols(); j++) {
      polys[i][j] = faces(i, j);
    }
  }

  SimplePolygonMesh pmesh(polys, coords);

  // Call the mesh writer
  pmesh.writeMesh(filename);
}

DenseMatrix<double> read_point_cloud(std::string filename) {

  std::unique_ptr<PointCloud> cloud;
  std::unique_ptr<PointPositionGeometry> geom;
  std::tie(cloud, geom) = readPointCloud(filename);

  return EigenMap<double, 3, Eigen::RowMajor>(geom->positions);
}

void write_point_cloud(DenseMatrix<double> points, std::string filename) {

  // Copy in to the point cloud object
  PointCloud cloud(points.rows());
  PointPositionGeometry geom(cloud);
  for (size_t i = 0; i < points.rows(); i++) {
    for (size_t j = 0; j < 3; j++) {
      geom.positions[i][j] = points(i, j);
    }
  }

  // Call the writer
  writePointCloud(cloud, geom, filename);
}


// Actual binding code
// clang-format off
void bind_io(py::module& m) {
  
  m.def("read_mesh", &read_mesh, "Read a mesh from file.", py::arg("filename"));
  m.def("write_mesh", &write_mesh, "Write a mesh to file.", py::arg("verts"), py::arg("faces"), py::arg("filename"));
  
  m.def("read_point_cloud", &read_point_cloud, "Read a point cloud from file.", py::arg("filename"));
  m.def("write_point_cloud", &write_point_cloud, "Write a point cloud to file.", py::arg("points"), py::arg("filename"));
}
