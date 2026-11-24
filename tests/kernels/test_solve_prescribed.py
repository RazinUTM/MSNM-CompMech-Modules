import numpy as np
import pytest

from core.solver.solve_prescribed import solve_prescribed


def test_multiple_prescribed_dofs_with_applied_support_loads():
    K = [[2, -2, 0], [-2, 5, -3], [0, -3, 3]]
    U, R = solve_prescribed(K, [2, 5, 1], [2, 0], [2, -1])
    np.testing.assert_allclose(U, [-1, 1.8, 2])
    np.testing.assert_allclose(R, [-0.4, -7.6])


def test_empty_partitions():
    np.testing.assert_allclose(solve_prescribed([[2]], [4], [], [])[0], [2])
    U, R = solve_prescribed([[2]], [4], [0], [3])
    np.testing.assert_allclose(U, [3])
    np.testing.assert_allclose(R, [2])


def test_duplicate_dof_is_rejected():
    with pytest.raises(ValueError, match="unique"):
        solve_prescribed([[2]], [4], [0, 0], [3, 3])
