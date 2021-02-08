#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/pointcloud/point_cloud.h"
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

/*
std::tuple<SparseMatrix<double>, SparseMatrix<double>>
buildMeshLaplacian(const DenseMatrix<double>& vMat, const DenseMatrix<size_t>& fMat, double mollifyFactor) {

  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = makeSurfaceMeshAndGeometry(vMat, fMat);

  SparseMatrix<double> L, M;
  std::tie(L, M) = buildTuftedLaplacian(*mesh, *geometry, mollifyFactor);

  return std::make_tuple(L, M);
}
*/

/*
std::tuple<SparseMatrix<double>, SparseMatrix<double>> buildPointCloudLaplacian(const DenseMatrix<double>& vMat,
                                                                                double mollifyFactor, size_t nNeigh) {

  SimplePolygonMesh cloudMesh;

  // Copy to std vector representation
  cloudMesh.vertexCoordinates.resize(vMat.rows());
  for (size_t iP = 0; iP < cloudMesh.vertexCoordinates.size(); iP++) {
    cloudMesh.vertexCoordinates[iP] = Vector3{vMat(iP, 0), vMat(iP, 1), vMat(iP, 2)};
  }

  // Generate the local triangulations for the point cloud
  Neighbors_t neigh = generate_knn(cloudMesh.vertexCoordinates, nNeigh);
  std::vector<Vector3> normals = generate_normals(cloudMesh.vertexCoordinates, neigh);
  std::vector<std::vector<Vector2>> coords = generate_coords_projection(cloudMesh.vertexCoordinates, normals, neigh);
  LocalTriangulationResult localTri = build_delaunay_triangulations(coords, neigh, false);

  // Take the union of all triangles in all the neighborhoods
  for (size_t iPt = 0; iPt < cloudMesh.vertexCoordinates.size(); iPt++) {
    const std::vector<size_t>& thisNeigh = neigh[iPt];
    size_t nNeigh = thisNeigh.size();

    // Accumulate over triangles
    for (const auto& tri : localTri.pointTriangles[iPt]) {
      std::array<size_t, 3> triGlobal = {thisNeigh[tri[0]], thisNeigh[tri[1]], thisNeigh[tri[2]]};
      cloudMesh.polygons.push_back({triGlobal[0], triGlobal[1], triGlobal[2]});
    }
  }


  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  std::tie(mesh, geometry) = makeSurfaceMeshAndGeometry(cloudMesh.polygons, cloudMesh.vertexCoordinates);

  SparseMatrix<double> L, M;
  std::tie(L, M) = buildTuftedLaplacian(*mesh, *geometry, mollifyFactor);

  L = L / 3.;
  M = M / 3.;

  return std::make_tuple(L, M);
}
*/

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


// Actual binding code
// clang-format off
PYBIND11_MODULE(potpourri3d_bindings, m) {
  m.doc() = "potpourri3d low-level bindings";
  
  m.def("read_mesh", &read_mesh, "Read a mesh from file.", py::arg("filename"));
  m.def("write_mesh", &write_mesh, "Write a mesh to file.", py::arg("verts"), py::arg("faces"), py::arg("filename"));
 
  /*
  m.def("buildPointCloudLaplacian", &buildPointCloudLaplacian, "build the point cloud Laplacian", 
      py::arg("vMat"), py::arg("mollifyFactor"), py::arg("nNeigh"));
      */
}

// clang-format on
