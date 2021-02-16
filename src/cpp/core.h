#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "Eigen/Dense"

namespace py = pybind11;

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace geometrycentral::pointcloud;


// Forward declare module builders in various source files
void bind_io(py::module& m);
void bind_mesh(py::module& m);
void bind_point_cloud(py::module& m);
