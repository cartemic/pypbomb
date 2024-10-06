import numpy as np

from pypbomb.data import _units, stress_limits
from pypbomb.data.stress_limits import PipingType


def test_max_stress():
    no_interp = stress_limits.max_stress(PipingType.seamless, "304", _units.f_to_k(500))
    with_interp = stress_limits.max_stress(PipingType.welded, "316", _units.f_to_k(550))
    assert np.isclose(no_interp, 84116060.8586595), "On temperature boundary"
    assert np.isclose(with_interp, 71360756.5481251), "Interpolation required"
