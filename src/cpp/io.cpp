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
#include <pybind11/stl.h>

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

std::tuple<DenseMatrix<double>, std::vector<std::vector<int64_t>>> read_polygon_mesh(std::string filename) {

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

  std::vector<std::vector<int64_t>> polygons(pmesh.nFaces());
  for (size_t i = 0; i < pmesh.nFaces(); i++) {
    size_t fDegree = pmesh.polygons[i].size();
    polygons[i].resize(fDegree);
    for (size_t j = 0; j < fDegree; j++) {
      polygons[i][j] = pmesh.polygons[i][j];
    }
  }

  return std::make_tuple(V, polygons);
}

namespace { // anonymous helers
SimplePolygonMesh buildMesh(const DenseMatrix<double>& verts, const DenseMatrix<int64_t>& faces,
                            const DenseMatrix<double>& corner_UVs) {
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
  std::vector<std::vector<Vector2>> corner_params;
  if (corner_UVs.size() > 0) {
    corner_params.resize(faces.rows());
    for (size_t i = 0; i < faces.rows(); i++) {
      corner_params[i].resize(faces.cols());
      for (size_t j = 0; j < faces.cols(); j++) {
        size_t ind = i * faces.cols() + j;
        for (size_t k = 0; k < 2; k++) {
          corner_params[i][j][k] = corner_UVs(ind, k);
        }
      }
    }
  }

  return SimplePolygonMesh(polys, coords, corner_params);
}
SimplePolygonMesh buildMesh(const DenseMatrix<double>& verts, const DenseMatrix<int64_t>& faces) {
  DenseMatrix<double> empty_UVs = DenseMatrix<double>::Zero(0, 2);
  return buildMesh(verts, faces, empty_UVs);
}
} // namespace

void write_mesh(DenseMatrix<double> verts, DenseMatrix<int64_t> faces, std::string filename) {
  SimplePolygonMesh pmesh = buildMesh(verts, faces);
  pmesh.writeMesh(filename);
}

void write_mesh_pervertex_uv(DenseMatrix<double> verts, DenseMatrix<int64_t> faces, DenseMatrix<double> UVs,
                             std::string filename) {
  size_t V = verts.rows();
  size_t F = faces.rows();
  size_t D = faces.cols();

  // expand out to per-corner UVs
  DenseMatrix<double> face_UVs = DenseMatrix<double>::Zero(F * D, 2);
  for (size_t i = 0; i < F; i++) {
    for (size_t j = 0; j < D; j++) {
      size_t vInd = faces(i, j);
      for (size_t k = 0; k < 2; k++) {
        face_UVs(i * D + j, k) = UVs(vInd, k);
      }
    }
  }

  SimplePolygonMesh pmesh = buildMesh(verts, faces, face_UVs);
  pmesh.writeMesh(filename);
}

void write_mesh_perface_uv(DenseMatrix<double> verts, DenseMatrix<int64_t> faces, DenseMatrix<double> UVs,
                           std::string filename) {

  size_t V = verts.rows();
  size_t F = faces.rows();
  size_t D = faces.cols();

  // expand out to per-corner UVs
  DenseMatrix<double> face_UVs = DenseMatrix<double>::Zero(F * D, 2);
  for (size_t i = 0; i < F; i++) {
    for (size_t j = 0; j < D; j++) {
      for (size_t k = 0; k < 2; k++) {
        face_UVs(i * D + j, k) = UVs(i, k);
      }
    }
  }

  SimplePolygonMesh pmesh = buildMesh(verts, faces, face_UVs);
  pmesh.writeMesh(filename);
}

void write_mesh_percorner_uv(DenseMatrix<double> verts, DenseMatrix<int64_t> faces, DenseMatrix<double> UVs,
                             std::string filename) {
  SimplePolygonMesh pmesh = buildMesh(verts, faces, UVs);
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
  m.def("read_polygon_mesh", &read_polygon_mesh, "Read a polygon mesh from file.", py::arg("filename"));

  m.def("write_mesh", &write_mesh, "Write a mesh to file.", py::arg("verts"), py::arg("faces"), py::arg("filename"));
  m.def("write_mesh_pervertex_uv", &write_mesh_pervertex_uv, "Write a mesh to file.", py::arg("verts"), py::arg("faces"), py::arg("UVs"), py::arg("filename"));
  m.def("write_mesh_perface_uv", &write_mesh_perface_uv, "Write a mesh to file.", py::arg("verts"), py::arg("faces"), py::arg("UVs"), py::arg("filename"));
  m.def("write_mesh_percorner_uv", &write_mesh_percorner_uv, "Write a mesh to file.", py::arg("verts"), py::arg("faces"), py::arg("UVs"), py::arg("filename"));
  
  m.def("read_point_cloud", &read_point_cloud, "Read a point cloud from file.", py::arg("filename"));
  m.def("write_point_cloud", &write_point_cloud, "Write a point cloud to file.", py::arg("points"), py::arg("filename"));
}
