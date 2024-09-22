import os

import numpy as np

from .. import thermochem
from ..thermochem import find_mechanisms


def test_laminar_flame_speed():
    initial_temperature = 300
    initial_pressure = 101325
    species = {
        "CH4": 0.095057034220532327,
        "O2": 0.19011406844106465,
        "N2": 0.71482889733840305,
    }
    # [1] C. K. Law, Combustion Physics. Cambridge: Cambridge University Press, 2006.
    expected = 0.39  # value approximated from [1] fig. 7.7.7
    test = thermochem.laminar_flame_speed(initial_temperature, initial_pressure, species, "gri30.yaml")

    assert np.isclose(expected, test, 0.05)


def test_sound_speed_eq():
    # check air at 1 atm and 20°C against ideal gas calculation
    gamma = 1.4
    rr = 8.31451
    tt = 293.15
    mm = 0.0289645
    ideal_gas_calc = np.sqrt(gamma * rr * tt / mm)

    temp = 293.15
    press = 101325
    species = {"O2": 1, "N2": 3.76}
    mechanism = "gri30.yaml"
    test = thermochem.sound_speed_eq(temp, press, species, mechanism)

    assert np.isclose(ideal_gas_calc, test, 0.005)


def test_reflected_cj_shock():
    # `reflected_cj_shock` integrates functions already tested as part of `sd`; this test just is to ensure
    # that the integration doesn't cause any errors
    initial_temperature = 300
    initial_pressure = 101325
    species_dict = {"H2": 1, "O2": 0.5}
    mechanism = "gri30.yaml"
    thermochem.reflected_cj_shock(
        initial_temperature=initial_temperature,
        initial_pressure=initial_pressure,
        species=species_dict,
        mechanism=mechanism,
        parallelize=True,
    )


def test_find_mechanisms():
    found = find_mechanisms()

    assert "gri30.yaml" in found.mechanisms
    assert os.path.exists(found.path)
