import numpy as np

from pypbomb.data import _units, flange
from pypbomb.data.flange import FlangePressureLimits


def test_all_available_materials():
    assert "304L" in flange.available_materials()


def test_available_materials_for_flange_class():
    assert "304L" in flange.available_materials(flange_class=150)


def test_available_classes():
    assert 2500 in flange.available_classes(material="316")


def test_flange_limits_max_pressure():
    fsl = FlangePressureLimits(
        temperature=np.array([0, 1]),  # temp is sorted ascending during db read
        pressure_limit=np.array([2, 1]),  # max pressure decreases with temp
    )
    assert np.isclose(fsl.max_pressure(0.5), 1.5)


def test_flange_limits_max_temperature():
    fsl = FlangePressureLimits(
        temperature=np.array([0, 1]),  # temp is sorted ascending during db read
        pressure_limit=np.array([2, 1]),  # max pressure decreases with temp
    )
    assert np.isclose(fsl.max_temperature(1.5), 0.5)


def test_pressure_limits():
    result = flange.pressure_limits("316L", 1500)
    assert np.allclose(
        result.temperature,
        _units.f_to_k(
            np.array(
                [
                    -29.0,
                    38.0,
                    50.0,
                    100.0,
                    150.0,
                    200.0,
                    250.0,
                    300.0,
                    325.0,
                    350.0,
                    375.0,
                    400.0,
                    425.0,
                    450.0,
                ]
            )
        ),
    )
    assert np.allclose(
        result.pressure_limit,
        _units.kpsi_to_pa(
            np.array(
                [
                    206.8,
                    206.8,
                    200.1,
                    173.9,
                    157,
                    145.8,
                    137.3,
                    130.3,
                    127.4,
                    125.4,
                    123.8,
                    121.5,
                    119.3,
                    117.1,
                ]
            )
        ),
    )


def test_get_class_good_output():
    temperature = 435.93  # 325 F -> K
    conditions = (
        # material, pressure, expected class
        (
            "304",
            625_285_538.9183401,
            900,
        ),  # 90.69 kpsi -> Pa (2.1) -- just below limit for 900
        (
            "304",
            625_423_434.0642036,
            1500,
        ),  # 90.71 kpsi -> Pa (2.1) -- just above limit for 900
    )

    for material, pressure, expected in conditions:
        test_class = flange.minimum_flange_class(temperature=temperature, pressure=pressure, material=material)

        assert test_class == expected
