import numpy as np
from cantera import Solution

from pypbomb import sd


def test_curve_fit():
    expected = sd.CurveFit(a=7.0, b=2.0, c=12.0, r2=1.0)
    x_values = np.linspace(0, 10, 100)
    y_values = expected.a * np.power(x_values, 2) + expected.b * x_values + expected.c
    test = sd.cj_curve_fit(x_values, y_values)

    assert np.isclose(test.a, expected.a)
    assert np.isclose(test.b, expected.b)
    assert np.isclose(test.c, expected.c)
    assert np.isclose(test.r2, expected.r2)


class TestCalculateCJSpeed:
    @staticmethod
    def test_serial():
        """
        Regression test against SDToolbox CJ speed calc (serial calculation)
        """
        expected = 2353.2706464533471
        test = sd.Detonation.cj_speed(
            initial_pressure=101325,
            initial_temperature=300,
            mole_fractions="H2:0.5333 O2:0.26667 AR:0.2",
            mechanism="gri30.yaml",
            parallelize=False,
            with_state=False,
        )

        assert np.isclose(expected, test.speed)
        assert test.state is None

    @staticmethod
    def test_parallel():
        """
        Regression test against SDToolbox CJ speed calc (parallel calculation)
        """
        expected = 2353.2706464533471
        test = sd.Detonation.cj_speed(
            initial_pressure=101325,
            initial_temperature=300,
            mole_fractions={"H2": 0.5333, "O2": 0.26667, "AR": 0.2},
            mechanism="gri30.yaml",
            parallelize=True,
            with_state=False,
        )

        assert np.isclose(expected, test.speed)
        assert test.state is None

    @staticmethod
    def test_with_state():
        test = sd.Detonation.cj_speed(
            initial_pressure=101325,
            initial_temperature=300,
            mole_fractions="H2:0.5333 O2:0.26667 AR:0.2",
            mechanism="gri30.yaml",
            parallelize=True,
            with_state=True,
        )

        assert isinstance(test.state, Solution)


def test_calculate_cj_state():
    # compare against SDToolbox results
    mechanism = "gri30.yaml"

    initial_gas = Solution(mechanism)
    initial_gas.TPX = 300, 101325, {"H2": 1}

    working_gas = Solution(mechanism)
    working_gas.TPX = 300, 101325 * 2, {"H2": 1}

    v_test = sd.Detonation.cj_state(working_gas, initial_gas, 1e-5, 1e-5, 1.5)
    v_expected = 1700.36
    # relaxed velocity requirement -- new version of cantera
    assert np.isclose(v_test, v_expected, 1e-4), "CJ velocity calculation out of spec"
    assert np.isclose(working_gas.T, 355.77590742266216, 1e-7), "CJ temperature calculation out of spec"
    assert np.isclose(working_gas.P, 180244.9690980063, 1e-7), "CJ pressure calculation out of spec"
    species_result = working_gas.mole_fraction_dict()
    for species, mf_expected in {"H": 2.8407416566652653e-30, "H2": 1.0}.items():
        err = f"{species} mole fraction calculation out of spec"
        assert np.isclose(species_result.get(species, np.nan), mf_expected, 1e-7), err
