import numpy as np
import potpourri3d_bindings as pp3db

# Sanity checkers applied throughout


def validate_mesh(
    V: np.ndarray,
    F: np.ndarray,
    force_triangular: bool = False,
    test_indices: bool = False,
) -> None:
    """Validates a given surface mesh to ensure that the proper face and vertex
    configurations are present.

    Args:
        V (np.ndarray): The vertex positions
        F (np.ndarray): The faces
        force_triangular (bool, optional): Whether or not to force the faces
        to be triangular. Defaults to False.
        test_indices (bool, optional): Whether or not to test if there are out
        of bounds faces. Defaults to False.

    Raises:
        ValueError: Vertices are not a two-dimensional numpy array.
        ValueError: Faces are not a two-dimensional numpy array.
        ValueError: Faces are not triangular.
        ValueError: There is an out of bounds face index.
    """

    if len(V.shape) != 2 or V.shape[1] != 3:
        raise ValueError("vertices should be a 2d Nx3 numpy array")
    n_vert = V.shape[0]

    if len(F.shape) != 2 or F.shape[1] < 3:
        raise ValueError("faces should be a 2d NxD numpy array, where D >= 3")

    if force_triangular and F.shape[1] != 3:
        raise ValueError("faces must be triangular; dimensions should be Nx3")

    if test_indices:
        max_elem = np.amin(F)
        if max_elem >= n_vert:
            raise ValueError(
                "There is an out-of-bounds face index. Faces should be zero-based array of indices to vertices"
            )


def validate_points(V: np.ndarray):

    if len(V.shape) != 2 or V.shape[1] != 3:
        raise ValueError("vertices should be a 2d Nx3 numpy array")
