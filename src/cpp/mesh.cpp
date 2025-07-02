#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/fast_marching_method.h"
#include "geometrycentral/surface/flip_geodesics.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/marching_triangles.h"
#include "geometrycentral/surface/mesh_graph_algorithms.h"
#include "geometrycentral/surface/polygon_mesh_heat_solver.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"

#include <memory>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Eigen/Dense"

#include "heat_helpers.h"

namespace py = pybind11;

using namespace geometrycentral;
using namespace geometrycentral::surface;


// For overloaded functions, with C++11 compiler only
template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

// Return the E x 2 integer-valued matrix representing the edges of a mesh as internal to geometry-central.
DenseMatrix<int64_t> edges(DenseMatrix<double> verts, DenseMatrix<int64_t> faces) {

  std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh(faces));
  DenseMatrix<int64_t> E(mesh->nEdges(), 2);
  for (size_t i = 0; i < mesh->nEdges(); i++) {
    E(i, 0) = mesh->edge(i).firstVertex().getIndex();
    E(i, 1) = mesh->edge(i).secondVertex().getIndex();
  }
  return E;
}

// Wrapper for FMMDistance that constructs an internal mesh and geometry
class FastMarchingDistanceEigen {

public:
  FastMarchingDistanceEigen(DenseMatrix<double> verts, DenseMatrix<int64_t> faces) {
    // Construct the internal mesh and geometry
    mesh.reset(new SurfaceMesh(faces));
    geom.reset(new VertexPositionGeometry(*mesh));
    for (size_t i = 0; i < mesh->nVertices(); i++) {
      for (size_t j = 0; j < 3; j++) {
        geom->inputVertexPositions[i][j] = verts(i, j);
      }
    }
    geom->requireVertexIndices();
    geom->requireEdgeIndices();
  }

  Vector<double> compute_distance(std::vector<std::vector<std::pair<int64_t, std::vector<double>>>> curves,
                                  std::vector<std::vector<double>> distances, bool sign = false) {

    size_t nCurves = curves.size();
    std::vector<std::vector<std::pair<SurfacePoint, double>>> initialDistances(nCurves);
    std::vector<std::vector<double>> initDists(nCurves);
    for (size_t i = 0; i < nCurves; i++) initDists[i] = std::vector<double>(curves[i].size(), 0.);

    for (size_t i = 0; i < std::min(distances.size(), nCurves); i++) {
      for (size_t j = 0; j < std::min(distances[i].size(), curves[i].size()); j++) {
        initDists[i][j] = distances[i][j];
      }
    }
    for (size_t i = 0; i < nCurves; i++) {
      for (size_t j = 0; j < curves[i].size(); j++) {
        initialDistances[i].emplace_back(toSurfacePoint(*mesh, curves[i][j]), initDists[i][j]);
      }
    }
    VertexData<double> phi = FMMDistance(*geom, initialDistances, sign);
    return phi.toVector();
  }

private:
  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
};

// Wrapper for marching triangles that constructs an internal mesh and geometry
class MarchingTrianglesEigen {

public:
  MarchingTrianglesEigen(DenseMatrix<double> verts, DenseMatrix<int64_t> faces) {
    // Construct the internal mesh and geometry
    mesh.reset(new SurfaceMesh(faces));
    geom.reset(new VertexPositionGeometry(*mesh));
    for (size_t i = 0; i < mesh->nVertices(); i++) {
      for (size_t j = 0; j < 3; j++) {
        geom->inputVertexPositions[i][j] = verts(i, j);
      }
    }

    geom->requireVertexIndices();
    geom->requireEdgeIndices();
    geom->requireFaceIndices();
  }

  std::vector<std::vector<std::pair<int64_t, std::vector<double>>>> marching_triangles(DenseMatrix<double> u,
                                                                                       double isoval) {

    VertexData<double> f(*mesh, u);
    std::vector<std::vector<SurfacePoint>> curves = marchingTriangles(*geom, f, isoval);
    size_t nCurves = curves.size();
    std::vector<std::vector<std::pair<int64_t, std::vector<double>>>> contour(nCurves);
    for (size_t i = 0; i < nCurves; i++) {
      std::pair<int64_t, std::vector<double>> pt;
      for (size_t j = 0; j < curves[i].size(); j++) {
        const SurfacePoint& p = curves[i][j];
        switch (p.type) {
        case (SurfacePointType::Vertex): {
          pt.first = geom->vertexIndices[p.vertex];
          break;
        }
        case (SurfacePointType::Edge): {
          pt.first = geom->edgeIndices[p.edge];
          pt.second = {p.tEdge};
          break;
        }
        default: {
          break;
        }
        }
        contour[i].push_back(pt);
      }
    }
    return contour;
  }

private:
  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
};

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
  VectorHeatMethodEigen(DenseMatrix<double> verts, DenseMatrix<int64_t> faces, double tCoef = 1.0,
                        bool useIntrinsicDelaunay_ = true)
      : useIntrinsicDelaunay(useIntrinsicDelaunay_) {

    // Construct the internal mesh and geometry
    mesh.reset(new ManifoldSurfaceMesh(faces));
    geom.reset(new VertexPositionGeometry(*mesh));
    for (size_t i = 0; i < mesh->nVertices(); i++) {
      for (size_t j = 0; j < 3; j++) {
        geom->inputVertexPositions[i][j] = verts(i, j);
      }
    }

    if (useIntrinsicDelaunay) {
      // Build the intrinsic triangulation, if using
      signpostTri.reset(new SignpostIntrinsicTriangulation(*mesh, *geom));
      signpostTri->flipToDelaunay();

      solver.reset(new VectorHeatMethodSolver(*signpostTri, tCoef));
    } else {
      solver.reset(new VectorHeatMethodSolver(*geom, tCoef));
    }
  }

  // Extend scalars from a collection of vertices
  Vector<double> extend_scalar(Vector<int64_t> sourceVerts, Vector<double> values) {
    std::vector<std::tuple<Vertex, double>> sources;
    for (size_t i = 0; i < sourceVerts.rows(); i++) {
      sources.emplace_back(get_compute_geom()->mesh.vertex(sourceVerts(i)), values(i));
      // ^^^ works on intrinsic if using, because vertices are unchanged
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
      // we always want to do this on the input mesh even if using intrinsic
      basisX[v] = geom->vertexTangentBasis[v][0];
      basisY[v] = geom->vertexTangentBasis[v][1];
    }

    return std::tuple<DenseMatrix<double>, DenseMatrix<double>, DenseMatrix<double>>(
        EigenMap<double, 3>(basisX), EigenMap<double, 3>(basisY), EigenMap<double, 3>(geom->vertexNormals));
  }

  SparseMatrix<std::complex<double>> get_connection_laplacian() {

    get_compute_geom()->requireVertexConnectionLaplacian();
    SparseMatrix<std::complex<double>> Lconn = get_compute_geom()->vertexConnectionLaplacian;
    get_compute_geom()->unrequireVertexConnectionLaplacian();
    return Lconn;
  }

  // TODO think about how to pass tangent frames around
  DenseMatrix<double> transport_tangent_vectors(Vector<int64_t> sourceVerts, DenseMatrix<double> values) {

    // Pack it as a Vector2
    std::vector<std::tuple<Vertex, Vector2>> sources;
    for (size_t i = 0; i < sourceVerts.rows(); i++) {
      sources.emplace_back(get_compute_geom()->mesh.vertex(sourceVerts(i)), Vector2{values(i, 0), values(i, 1)});
      // ^^^ works on intrinsic if using, because vertices and their tangent spaces are unchanged
    }
    VertexData<Vector2> ext = solver->transportTangentVectors(sources);

    return EigenMap<double, 2>(ext);
  }

  DenseMatrix<double> transport_tangent_vector(int64_t sourceVert, DenseMatrix<double> values) {

    // Pack it as a Vector2
    std::vector<std::tuple<Vertex, Vector2>> sources;
    sources.emplace_back(get_compute_geom()->mesh.vertex(sourceVert), Vector2{values(0), values(1)});
    VertexData<Vector2> ext = solver->transportTangentVectors(sources);

    return EigenMap<double, 2>(ext);
  }


  DenseMatrix<double> compute_log_map(int64_t sourceVert, std::string strategy) {
    LogMapStrategy strategyEnum = toLogmapStrategy(strategy);
    return EigenMap<double, 2>(solver->computeLogMap(get_compute_geom()->mesh.vertex(sourceVert), strategyEnum));
  }

private:
  const bool useIntrinsicDelaunay;
  std::unique_ptr<SignpostIntrinsicTriangulation> signpostTri;

  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::unique_ptr<VectorHeatMethodSolver> solver;

  // helpers
  IntrinsicGeometryInterface* get_compute_geom() {
    return useIntrinsicDelaunay ? static_cast<IntrinsicGeometryInterface*>(signpostTri.get())
                                : static_cast<IntrinsicGeometryInterface*>(geom.get());
  }
};

// A wrapper class for the signed heat method solver, which exposes Eigen in/out
class SignedHeatMethodEigen {

  // TODO use intrinsic triangulations here

public:
  SignedHeatMethodEigen(DenseMatrix<double> verts, DenseMatrix<int64_t> faces, double tCoef = 1.0) {

    // Construct the internal mesh and geometry
    mesh.reset(new SurfaceMesh(faces));
    geom.reset(new VertexPositionGeometry(*mesh));
    for (size_t i = 0; i < mesh->nVertices(); i++) {
      for (size_t j = 0; j < 3; j++) {
        geom->inputVertexPositions[i][j] = verts(i, j);
      }
    }

    // Build the solver
    solver.reset(new SignedHeatSolver(*geom, tCoef));
  }

  Vector<double> compute_distance(const std::vector<std::vector<std::pair<int64_t, std::vector<double>>>>& curves,
                                  const std::vector<bool>& isSigned,
                                  const std::vector<std::pair<int64_t, std::vector<double>>>& points,
                                  bool preserveSourceNormals = false, const std::string& levelSetConstraint = "ZeroSet",
                                  double softLevelSetWeight = -1.0) {

    std::vector<Curve> sourceCurves = toSignedCurves(*mesh, curves, isSigned);
    std::vector<SurfacePoint> sourcePoints = toSurfacePoints(*mesh, points);
    SignedHeatOptions solveOptions = toSignedHeatOptions(preserveSourceNormals, levelSetConstraint, softLevelSetWeight);
    VertexData<double> phi = solver->computeDistance(sourceCurves, sourcePoints, solveOptions);
    return phi.toVector();
  }

private:
  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::unique_ptr<SignedHeatSolver> solver;
};

// A wrapper class for the polygon mesh heat method solver, which exposes Eigen in/out
class PolygonMeshHeatMethodEigen {

public:
  PolygonMeshHeatMethodEigen(const DenseMatrix<double>& verts, const std::vector<std::vector<size_t>>& faces,
                             double tCoef = 1.0) {

    // Construct the internal mesh and geometry
    mesh.reset(new SurfaceMesh(faces));
    geom.reset(new VertexPositionGeometry(*mesh));
    for (size_t i = 0; i < mesh->nVertices(); i++) {
      for (size_t j = 0; j < 3; j++) {
        geom->inputVertexPositions[i][j] = verts(i, j);
      }
    }

    // Build the solver
    solver.reset(new PolygonMeshHeatSolver(*geom, tCoef));
  }

  Vector<double> compute_distance(Vector<int64_t> sourceVerts) {
    std::vector<Vertex> sources;
    for (size_t i = 0; i < sourceVerts.rows(); i++) {
      sources.push_back(mesh->vertex(sourceVerts(i)));
    }
    VertexData<double> dist = solver->computeDistance(sources);
    return dist.toVector();
  }

  Vector<double> extend_scalar(Vector<int64_t> sourceVerts, Vector<double> values) {
    std::vector<std::tuple<Vertex, double>> sources;
    for (size_t i = 0; i < sourceVerts.rows(); i++) {
      sources.emplace_back(mesh->vertex(sourceVerts(i)), values(i));
    }
    VertexData<double> ext = solver->extendScalars(sources);
    return ext.toVector();
  }

  // Returns an extrinsic representation of the tangent frame being used internally, as X/Y/N vectors.
  // This is a direct copy of what's already in the VectorHeatMethodEigen class. It's not ideal to repeat code, but this
  // also seems like the simplest solution without making breaking changes.
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

  DenseMatrix<double> transport_tangent_vectors(Vector<int64_t> sourceVerts, DenseMatrix<double> values) {
    std::vector<std::tuple<Vertex, Vector2>> sources;
    for (size_t i = 0; i < sourceVerts.rows(); i++) {
      sources.emplace_back(mesh->vertex(sourceVerts(i)), Vector2{values(i, 0), values(i, 1)});
    }
    VertexData<Vector2> ext = solver->transportTangentVectors(sources);

    return EigenMap<double, 2>(ext);
  }

  Vector<double> compute_signed_distance(const std::vector<std::vector<int64_t>>& curves,
                                         const std::string& levelSetConstraint = "ZeroSet") {

    std::vector<std::vector<Vertex>> sourceCurves = toSignedVertices(*mesh, curves);
    SignedHeatOptions solveOptions = toSignedHeatOptions(false, levelSetConstraint, -1.0);
    VertexData<double> dist = solver->computeSignedDistance(sourceCurves, solveOptions.levelSetConstraint);
    return dist.toVector();
  }


private:
  std::unique_ptr<SurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::unique_ptr<PolygonMeshHeatSolver> solver;
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
  DenseMatrix<double> find_geodesic_path(int64_t startVert, int64_t endVert, size_t maxIterations = INVALID_IND,
                                         double maxRelativeLengthDecrease = 0.) {

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
    flipNetwork->iterativeShorten(maxIterations, maxRelativeLengthDecrease);

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
  DenseMatrix<double> find_geodesic_path_poly(std::vector<int64_t> verts, size_t maxIterations = INVALID_IND,
                                              double maxRelativeLengthDecrease = 0.) {

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
    flipNetwork->iterativeShorten(maxIterations, maxRelativeLengthDecrease);

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
  DenseMatrix<double> find_geodesic_loop(std::vector<int64_t> verts, size_t maxIterations = INVALID_IND,
                                         double maxRelativeLengthDecrease = 0.) {

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
    flipNetwork->iterativeShorten(maxIterations, maxRelativeLengthDecrease);

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

private:
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
  std::unique_ptr<FlipEdgeNetwork> flipNetwork;
};

class GeodesicTracer {

public:
  GeodesicTracer(DenseMatrix<double> verts, DenseMatrix<int64_t> faces) {

    // Construct the internal mesh and geometry
    mesh.reset(new ManifoldSurfaceMesh(faces));
    geom.reset(new VertexPositionGeometry(*mesh));
    for (size_t i = 0; i < mesh->nVertices(); i++) {
      for (size_t j = 0; j < 3; j++) {
        geom->inputVertexPositions[i][j] = verts(i, j);
      }
    }

    geom->requireVertexTangentBasis();
    geom->requireFaceTangentBasis();
  }

  // Generate a geodesic by tracing from a vertex along a tangent direction
  DenseMatrix<double> trace_geodesic_worker(SurfacePoint start_point, Vector2 start_dir,
                                            size_t max_iters = INVALID_IND) {

    TraceOptions opts;
    opts.includePath = true;
    opts.errorOnProblem = false;
    opts.barrierEdges = nullptr;
    opts.maxIters = max_iters;

    TraceGeodesicResult result = traceGeodesic(*geom, start_point, start_dir, opts);

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

  // Generate a geodesic by tracing from a vertex along a tangent direction
  DenseMatrix<double> trace_geodesic_from_vertex(int64_t startVert, Eigen::Vector3d direction_xyz,
                                                 size_t max_iters = INVALID_IND) {

    // Project the input direction onto the tangent basis
    Vertex vert = mesh->vertex(startVert);
    Vector3 direction{direction_xyz(0), direction_xyz(1), direction_xyz(2)};
    Vector3 bX = geom->vertexTangentBasis[vert][0];
    Vector3 bY = geom->vertexTangentBasis[vert][1];

    // magnitude matters! it determines the length
    Vector2 tangent_dir{dot(direction, bX), dot(direction, bY)};

    return trace_geodesic_worker(SurfacePoint(vert), tangent_dir, max_iters);
  }

  // Generate a geodesic by tracing from a face along a tangent direction
  DenseMatrix<double> trace_geodesic_from_face(int64_t startFace, Eigen::Vector3d bary_coords,
                                               Eigen::Vector3d direction_xyz, size_t max_iters = INVALID_IND) {

    // Project the input direction onto the tangent basis
    Face face = mesh->face(startFace);
    Vector3 bary{bary_coords(0), bary_coords(1), bary_coords(2)};
    Vector3 direction{direction_xyz(0), direction_xyz(1), direction_xyz(2)};
    Vector3 bX = geom->faceTangentBasis[face][0];
    Vector3 bY = geom->faceTangentBasis[face][1];

    // magnitude matters! it determines the length
    Vector2 tangent_dir{dot(direction, bX), dot(direction, bY)};

    return trace_geodesic_worker(SurfacePoint(face, bary), tangent_dir, max_iters);
  }

private:
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geom;
};


// Actual binding code
// clang-format off
void bind_mesh(py::module& m) {

  m.def("edges", &edges, py::arg("verts"), py::arg("faces"));

  py::class_<FastMarchingDistanceEigen>(m, "MeshFastMarchingDistance")
        .def(py::init<DenseMatrix<double>, DenseMatrix<int64_t>>())
        .def("compute_distance", &FastMarchingDistanceEigen::compute_distance, py::arg("curves"), py::arg("distances"), py::arg("sign"));

  py::class_<MarchingTrianglesEigen>(m, "MeshMarchingTriangles")
        .def(py::init<DenseMatrix<double>, DenseMatrix<int64_t>>())
        .def("marching_triangles", &MarchingTrianglesEigen::marching_triangles, py::arg("u"), py::arg("isoval"));

  py::class_<HeatMethodDistanceEigen>(m, "MeshHeatMethodDistance")
        .def(py::init<DenseMatrix<double>, DenseMatrix<int64_t>, double, bool>())
        .def("compute_distance", &HeatMethodDistanceEigen::compute_distance, py::arg("source_vert"))
        .def("compute_distance_multisource", &HeatMethodDistanceEigen::compute_distance_multisource, py::arg("source_verts"));
 

  py::class_<VectorHeatMethodEigen>(m, "MeshVectorHeatMethod")
        .def(py::init<DenseMatrix<double>, DenseMatrix<int64_t>, double, bool>())
        .def("extend_scalar", &VectorHeatMethodEigen::extend_scalar, py::arg("source_verts"), py::arg("values"))
        .def("get_tangent_frames", &VectorHeatMethodEigen::get_tangent_frames)
        .def("get_connection_laplacian", &VectorHeatMethodEigen::get_connection_laplacian)
        .def("transport_tangent_vector", &VectorHeatMethodEigen::transport_tangent_vector, py::arg("source_vert"), py::arg("vector"))
        .def("transport_tangent_vectors", &VectorHeatMethodEigen::transport_tangent_vectors, py::arg("source_verts"), py::arg("vectors"))
        .def("compute_log_map", &VectorHeatMethodEigen::compute_log_map, py::arg("source_vert"), py::arg("strategy"));

  py::class_<SignedHeatMethodEigen>(m, "MeshSignedHeatMethod")
        .def(py::init<DenseMatrix<double>, DenseMatrix<int64_t>, double>())
        .def("compute_distance", &SignedHeatMethodEigen::compute_distance,
             py::arg("curves") = std::vector<std::vector<std::pair<int64_t, std::vector<double>>>>(),
             py::arg("is_signed") = std::vector<bool>(),
             py::arg("points") = std::vector<std::pair<int64_t, std::vector<double>>>(),
             py::arg("preserve_source_normals") = false,
             py::arg("level_set_constraint") = "ZeroSet",
             py::arg("soft_level_set_weight") = -1.0);

  py::class_<PolygonMeshHeatMethodEigen>(m, "PolygonMeshHeatSolver")
        .def(py::init<DenseMatrix<double>, std::vector<std::vector<size_t>>, double>())
        .def("compute_distance", &PolygonMeshHeatMethodEigen::compute_distance, py::arg("source_verts"))
        .def("extend_scalar", &PolygonMeshHeatMethodEigen::extend_scalar, py::arg("source_verts"), py::arg("values"))
        .def("get_tangent_frames", &PolygonMeshHeatMethodEigen::get_tangent_frames)
        .def("transport_tangent_vectors", &PolygonMeshHeatMethodEigen::transport_tangent_vectors, py::arg("source_verts"), py::arg("vectors"))
        .def("compute_signed_distance", &PolygonMeshHeatMethodEigen::compute_signed_distance, 
             py::arg("curves") = std::vector<std::vector<int64_t>>(),
             py::arg("level_set_constraint") = "ZeroSet");

  py::class_<EdgeFlipGeodesicsManager>(m, "EdgeFlipGeodesicsManager")
        .def(py::init<DenseMatrix<double>, DenseMatrix<int64_t>>())
        .def("find_geodesic_path", &EdgeFlipGeodesicsManager::find_geodesic_path, py::arg("source_vert"), py::arg("target_vert"), py::arg("maxIterations"), py::arg("maxRelativeLengthDecrease"))
        .def("find_geodesic_path_poly", &EdgeFlipGeodesicsManager::find_geodesic_path_poly, py::arg("vert_list"), py::arg("maxIterations"), py::arg("maxRelativeLengthDecrease"))
        .def("find_geodesic_loop", &EdgeFlipGeodesicsManager::find_geodesic_loop, py::arg("vert_list"), py::arg("maxIterations"), py::arg("maxRelativeLengthDecrease"));


  py::class_<GeodesicTracer>(m, "GeodesicTracer")
        .def(py::init<DenseMatrix<double>, DenseMatrix<int64_t>>())
        .def("trace_geodesic_from_vertex", &GeodesicTracer::trace_geodesic_from_vertex, py::arg("start_vert"), py::arg("direction_xyz"), py::arg("max_iters"))
        .def("trace_geodesic_from_face", &GeodesicTracer::trace_geodesic_from_face, py::arg("start_face"), py::arg("bary_coords"), py::arg("direction_xyz"), py::arg("max_iters"));

}
