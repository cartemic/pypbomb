# -*- coding: utf-8 -*-

import numpy as np
import pytest
from cantera import Solution

from .. import sd


def test_curve_fit():
    abc = np.array([7.0, 2.0, 12.0, 1.0])

    x_values = np.linspace(0, 10, 100)
    y_values = abc[0] * np.power(x_values, 2) + abc[1] * x_values + abc[2]

    fit = sd.cj_curve_fit(x_values, y_values)
    assert np.allclose(abc, fit)


class TestCalculateCJSpeed:
    @staticmethod
    def test_no_multiprocessing():
        # Check for correct calculation when multiprocessing is not used
        test_args = [
            101325,
            300,
            "H2:0.5333 O2:0.26667 AR:0.2",
            "gri30.cti",
            False,  # use multiprocessing
            False,  # return r squared
            False,  # return state
        ]

        # CJ speed calculated by SDToolbox given *test_args
        original_cj_speed = 2353.2706464533471

        test_speed = sd.Detonation.cj_speed(*test_args)["cj speed"]

        assert np.allclose(original_cj_speed, test_speed)

    @staticmethod
    def test_multiprocessing():
        # Check for correct calculation when multiprocessing is used
        test_args = [
            101325,
            300,
            "H2:0.5333 O2:0.26667 AR:0.2",
            "gri30.cti",
            True,  # use multiprocessing
            False,  # return r squared
            False,  # return state
        ]

        # CJ speed calculated by SDToolbox given *test_args
        original_cj_speed = 2353.2706464533471

        test_speed = sd.Detonation.cj_speed(*test_args)["cj speed"]

        assert np.allclose(original_cj_speed, test_speed)

    @staticmethod
    def test_return_state():
        test_args = [
            101325,
            300,
            "H2:0.5333 O2:0.26667 AR:0.2",
            "gri30.cti",
            True,  # use multiprocessing
            False,  # return r squared
            True,  # return state
        ]
        test_info = sd.Detonation.cj_speed(*test_args)
        required_keys = {"cj speed", "cj state"}

        check_keys = len(required_keys.difference(test_info.keys())) == 0
        check_state = isinstance(test_info["cj state"], Solution)

        assert all([check_keys, check_state])

    @staticmethod
    def test_return_r_squared():
        test_args = [
            101325,
            300,
            "H2:0.5333 O2:0.26667 AR:0.2",
            "gri30.cti",
            True,  # use multiprocessing
            True,  # return r squared
            False,  # return state
        ]
        test_info = sd.Detonation.cj_speed(*test_args)
        required_keys = {"cj speed", "R^2"}

        check_keys = len(required_keys.difference(test_info.keys())) == 0
        check_r_squared = test_info["R^2"] >= 0.99

        assert all([check_keys, check_r_squared])

    @staticmethod
    def test_return_state_and_r_squared():
        test_args = [
            101325,
            300,
            "H2:0.5333 O2:0.26667 AR:0.2",
            "gri30.cti",
            True,  # use multiprocessing
            True,  # return r squared
            True,  # return state
        ]
        test_info = sd.Detonation.cj_speed(*test_args)
        required_keys = {"cj speed", "cj state", "R^2"}

        check_keys = len(required_keys.difference(test_info.keys())) == 0
        check_state = isinstance(test_info["cj state"], Solution)
        check_r_squared = test_info["R^2"] >= 0.99

        assert all([check_keys, check_state, check_r_squared])


class TestCalculateCJState:
    @staticmethod
    def test_good_input():
        # compare against SDToolbox results
        mechanism = "gri30.cti"

        initial_gas = Solution(mechanism)
        initial_gas.TPX = 300, 101325, {"H2": 1}

        working_gas = Solution(mechanism)
        working_gas.TPX = 300, 101325 * 2, {"H2": 1}

        cj_calcs = sd.Detonation.cj_state(working_gas, initial_gas, 1e-5, 1e-5, 1.5)

        good_temp = 355.77590742266216
        good_press = 180244.9690980063
        good_species = {"H": 2.8407416566652653e-30, "H2": 1.0}

        test_temp = cj_calcs[0].T
        check_temp = abs(test_temp - good_temp) / good_temp < 1e-7

        test_press = cj_calcs[0].P
        check_press = abs(test_press - good_press) / good_press < 1e-7

        # relaxed velocity requirement -- new version of cantera
        good_velocity = 1700.36
        test_velocity = cj_calcs[1]
        check_velocity = abs(test_velocity - good_velocity) / good_velocity <= 1e-1

        checks = [check_temp, check_press, check_velocity]

        # make sure the species in each solution are the same
        test_species = cj_calcs[0].mole_fraction_dict()
        for species in good_species:
            checks.append(good_species[species] - test_species[species] < 1e-7)

        assert all(checks)

    @staticmethod
    def test_no_convergence():
        # ensure the proper warning is generated when solution doesn't converge
        mechanism = "gri30.cti"

        initial_gas = Solution(mechanism)
        initial_gas.TPX = 300, 101325, {"H2": 1}

        working_gas = Solution(mechanism)
        working_gas.TPX = 300, 101325 * 2, {"H2": 1}

        with pytest.warns(Warning, match="No convergence within 1 iterations"):
            sd.Detonation.cj_state(working_gas, initial_gas, 1e-50, 1e-50, 1.5, 1)
