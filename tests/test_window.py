import numpy as np

from pypbomb import window


def test_safety_factor():
    expected = 4
    result = window.safety_factor(
        length=0.02,  # 20 mm -> m
        width=0.05,  # 50 mm -> m
        thickness=0.0012,  # 1.2 mm -> m
        internal_pressure=101_325,
        external_pressure=0,
        rupture_modulus=36542223.1599095,  # 5300 psi -> Pa
    )

    assert np.isclose(result, expected, rtol=0.01)


def test_minimum_thickness():
    expected = 0.0012  # 1.2 mm -> m
    result = window.minimum_thickness(
        length=0.02,  # 20 mm -> m
        width=0.05,  # 50 mm -> m
        fs=4,
        internal_pressure=101_325,
        external_pressure=0,
        rupture_modulus=36542223.1599095,  # 5300 psi -> Pa
    )
    assert np.isclose(result, expected, rtol=0.01)
