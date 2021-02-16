#include "core.h"


// Actual binding code
PYBIND11_MODULE(potpourri3d_bindings, m) {
  m.doc() = "potpourri3d low-level bindings";

  bind_io(m);
  bind_mesh(m);
  bind_point_cloud(m);
}
