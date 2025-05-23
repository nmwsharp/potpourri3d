#include "heat_helpers.h"

std::vector<Curve> toSignedCurves(SurfaceMesh& mesh,
                                  const std::vector<std::vector<std::pair<int64_t, std::vector<double>>>>& pythonCurves,
                                  const std::vector<bool>& isSigned) {

  std::vector<Curve> curves;

  for (size_t i = 0; i < pythonCurves.size(); i++) {
    Curve curve;
    curve.isSigned = (i < isSigned.size()) ? isSigned[i] : true; // default to signed if not specified

    for (const auto& node : pythonCurves[i]) {
      size_t elementIndex = node.first;
      const std::vector<double>& baryCoords = node.second;

      if (baryCoords.size() < 1) {
        curve.nodes.emplace_back(mesh.vertex(elementIndex));
      } else if (baryCoords.size() == 1) {
        curve.nodes.emplace_back(mesh.edge(elementIndex), baryCoords[0]);
      } else if (baryCoords.size() >= 2) {
        double tC = (baryCoords.size() == 3) ? baryCoords[2] : 1.0 - baryCoords[0] - baryCoords[1];
        Vector3 faceCoords = {baryCoords[0], baryCoords[1], tC};
        curve.nodes.emplace_back(mesh.face(elementIndex), faceCoords);
      } else {
        throw std::runtime_error("Invalid barycentric coordinates for a curve node.");
      }
    }
    curves.push_back(curve);
  }

  return curves;
}

std::vector<Curve> toSignedCurves(SurfaceMesh& mesh, const std::vector<std::vector<int64_t>>& pythonCurves,
                                  const std::vector<bool>& isSigned) {

  std::vector<Curve> curves;

  for (size_t i = 0; i < pythonCurves.size(); i++) {
    Curve curve;
    curve.isSigned = (i < isSigned.size()) ? isSigned[i] : true; // default to signed if not specified
    curve.nodes = toSurfacePoints(mesh, pythonCurves[i]);
    curves.push_back(curve);
  }

  return curves;
}

std::vector<std::vector<Vertex>> toSignedVertices(SurfaceMesh& mesh,
                                                  const std::vector<std::vector<int64_t>>& pythonCurves) {
  std::vector<std::vector<Vertex>> curves;
  for (size_t i = 0; i < pythonCurves.size(); i++) {
    curves.emplace_back();
    for (const auto& vIdx : pythonCurves[i]) {
      curves.back().push_back(mesh.vertex(vIdx));
    }
  }
  return curves;
}

std::vector<SurfacePoint> toSurfacePoints(SurfaceMesh& mesh,
                                          const std::vector<std::pair<int64_t, std::vector<double>>>& pythonPoints) {

  std::vector<SurfacePoint> points;

  for (const auto& point : pythonPoints) {
    size_t elementIndex = point.first;
    const std::vector<double>& baryCoords = point.second;

    if (baryCoords.size() < 1) {
      points.emplace_back(mesh.vertex(elementIndex));
    } else if (baryCoords.size() == 1) {
      points.emplace_back(mesh.edge(elementIndex), baryCoords[0]);
    } else if (baryCoords.size() >= 2) {
      double tC = (baryCoords.size() == 3) ? baryCoords[2] : 1.0 - baryCoords[0] - baryCoords[1];
      Vector3 faceCoords = {baryCoords[0], baryCoords[1], tC};
      points.emplace_back(mesh.face(elementIndex), faceCoords);
    } else {
      throw std::runtime_error("Invalid barycentric coordinates for surface point.");
    }
  }

  return points;
}

std::vector<SurfacePoint> toSurfacePoints(SurfaceMesh& mesh, const std::vector<int64_t>& pythonPoints) {

  std::vector<SurfacePoint> points;

  for (const auto& vIdx : pythonPoints) {
    points.emplace_back(mesh.vertex(vIdx));
  }

  return points;
}

SignedHeatOptions toSignedHeatOptions(bool preserveSourceNormals, const std::string& levelSetConstraint,
                                      double softLevelSetWeight) {

  auto toLower = [&](const std::string& s) -> std::string {
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return std::tolower(c); });
    return t;
  };

  SignedHeatOptions options;
  options.preserveSourceNormals = preserveSourceNormals;
  if (toLower(levelSetConstraint) == "none") {
    options.levelSetConstraint = LevelSetConstraint::None;
  }
  if (toLower(levelSetConstraint) == "zeroset") {
    options.levelSetConstraint = LevelSetConstraint::ZeroSet;
  }
  if (toLower(levelSetConstraint) == "multiple") {
    options.levelSetConstraint = LevelSetConstraint::Multiple;
  }
  options.softLevelSetWeight = softLevelSetWeight;
  return options;
}