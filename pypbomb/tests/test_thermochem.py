# -*- coding: utf-8 -*-

import os

import numpy as np
import pint

from .. import thermochem

_U = pint.UnitRegistry()
_Q = _U.Quantity


def test_calculate_laminar_flame_speed():
    initial_temperature = _Q(300, "K")
    initial_pressure = _Q(1, "atm")
    species = {
        "CH4": 0.095057034220532327,
        "O2": 0.19011406844106465,
        "N2": 0.71482889733840305,
    }
    good_result = 0.39  # value approximated from Law fig. 7.7.7
    test_flame_speed = thermochem.calculate_laminar_flame_speed(
        initial_temperature, initial_pressure, species, "gri30.cti"
    )
    assert abs(test_flame_speed.magnitude - good_result) / good_result < 0.05


def test_get_eq_sound_speed():
    # check air at 1 atm and 20°C against ideal gas calculation
    gamma = 1.4
    rr = 8.31451
    tt = 293.15
    mm = 0.0289645
    c_ideal = np.sqrt(gamma * rr * tt / mm)

    temp = _Q(20, "degC")
    press = _Q(1, "atm")
    species = {"O2": 1, "N2": 3.76}
    mechanism = "gri30.cti"
    c_test = thermochem.get_eq_sound_speed(temp, press, species, mechanism)

    assert abs(c_ideal - c_test.to("m/s").magnitude) / c_ideal <= 0.005


def test_calculate_reflected_shock_state():
    # this is just a handler for some sd2 functions, and this test is to ensure
    # that it doesn't throw any errors
    initial_temperature = _Q(80, "degF")
    initial_pressure = _Q(1, "atm")
    species_dict = {"H2": 1, "O2": 0.5}
    mechanism = "gri30.cti"
    thermochem.calculate_reflected_shock_state(initial_temperature, initial_pressure, species_dict, mechanism, _U)


class TestFindMechanisms:
    def test_mechs_only(self):
        assert "gri30.cti" in thermochem.find_mechanisms()

    def test_return_directory(self):
        checks = [False, False]
        mechs, path = thermochem.find_mechanisms(True)
        checks[0] = "gri30.cti" in mechs
        checks[1] = os.path.exists(path)
        assert all(checks)
