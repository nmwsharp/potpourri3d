#include "geometrycentral/surface/signed_heat_method.h"
#include "geometrycentral/surface/vector_heat_method.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// Convert from intermediary representation -- used for easier Python bindings -- to the input structs used in
// geometry-central.

SurfacePoint toSurfacePoint(SurfaceMesh& mesh, const std::pair<int64_t, std::vector<double>>& p);

std::vector<Curve> toSignedCurves(SurfaceMesh& mesh,
                                  const std::vector<std::vector<std::pair<int64_t, std::vector<double>>>>& pythonCurves,
                                  const std::vector<bool>& isSigned);

std::vector<Curve> toSignedCurves(SurfaceMesh& mesh, const std::vector<std::vector<int64_t>>& pythonCurves,
                                  const std::vector<bool>& isSigned);

std::vector<std::vector<Vertex>> toSignedVertices(SurfaceMesh& mesh,
                                                  const std::vector<std::vector<int64_t>>& pythonCurves);

std::vector<SurfacePoint> toSurfacePoints(SurfaceMesh& mesh,
                                          const std::vector<std::pair<int64_t, std::vector<double>>>& pythonPoints);

std::vector<SurfacePoint> toSurfacePoints(SurfaceMesh& mesh, const std::vector<int64_t>& pythonPoints);

SignedHeatOptions toSignedHeatOptions(bool preserveSourceNormals = false,
                                      const std::string& levelSetConstraint = "ZeroSet",
                                      double softLevelSetWeight = -1.);


LogMapStrategy toLogmapStrategy(std::string strategyName);