# -*- coding: utf-8 -*-
"""
This module contains functions for performing thermochemical calculations using
``cantera`` and ``pypbomb.sd``.
"""

import os

import cantera as ct
import numpy as np
import pint

from . import units, sd


_U = pint.UnitRegistry()


def calculate_laminar_flame_speed(
        initial_temperature,
        initial_pressure,
        species,
        mechanism,
        phase_specification="",
        unit_registry=_U
):
    """
    This function uses cantera to calculate the laminar flame speed of a given
    gas mixture.

    Parameters
    ----------
    initial_temperature : pint.Quantity
        Mixture initial temperature
    initial_pressure : pint.Quantity
        Mixture initial pressure
    species : dict or str
        Species definition for cantera
    mechanism : str
        String of mechanism to use (e.g. ``gri30.cti``)
    phase_specification : str, optional
        Phase specification for cantera solution
    unit_registry : pint.UnitRegistry, optional
        Unit registry for managing units to prevent conflicts with parent
        unit registry

    Returns
    -------
    pint.Quantity
        Laminar flame speed
    """
    gas = ct.Solution(mechanism, phase_specification)
    quant = unit_registry.Quantity

    initial_pressure = units.parse_quant_input(
        initial_pressure,
        unit_registry
    )
    initial_temperature = units.parse_quant_input(
        initial_temperature,
        unit_registry
    )
    units.check_pint_quantity(
        initial_pressure,
        "pressure",
        ensure_positive=True
    )
    units.check_pint_quantity(
        initial_temperature,
        "temperature",
        ensure_positive=True
    )

    gas.TPX = (
        initial_temperature.to("K").magnitude,
        initial_pressure.to("Pa").magnitude,
        species
    )

    # find laminar flame speed
    flame = ct.FreeFlame(gas)
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
    flame.solve(loglevel=0)

    return quant(flame.u[0], "m/s")


# noinspection SpellCheckingInspection
def get_eq_sound_speed(
        temperature,
        pressure,
        species,
        mechanism,
        phase_specification="",
        unit_registry=_U
):
    """
    Calculates the equilibrium speed of sound in a mixture

    Parameters
    ----------
    temperature : pint.Quantity or Tuple[float, str]
        Mixture initial temperature
    pressure : pint.Quantity or Tuple[float, str]
        Mixture initial pressure
    species : dict or str
        Species definition for cantera
    mechanism : str
        Desired chemical mechanism
    phase_specification : str, optional
        Phase specification for cantera solution
    unit_registry : pint.UnitRegistry, optional
        Unit registry for managing units to prevent conflicts with parent
        unit registry

    Returns
    -------
    sound_speed : pint.Quantity or Tuple[float, str]
        local speed of sound
    """
    quant = unit_registry.Quantity

    pressure = units.parse_quant_input(
        pressure,
        unit_registry
    )
    temperature = units.parse_quant_input(
        temperature,
        unit_registry
    )
    units.check_pint_quantity(
        pressure,
        "pressure",
        ensure_positive=True
    )
    units.check_pint_quantity(
        temperature,
        "temperature",
        ensure_positive=True
    )

    working_gas = ct.Solution(mechanism, phase_specification)
    working_gas.TPX = [
        temperature.to("K").magnitude,
        pressure.to("Pa").magnitude,
        species
        ]

    pressures = np.zeros(2)
    densities = np.zeros(2)

    # equilibrate gas at input conditions and collect pressure, density
    working_gas.equilibrate("TP")
    pressures[0] = working_gas.P
    densities[0] = working_gas.density

    # perturb pressure and equilibrate with constant P, s to get dp/drho|s
    pressures[1] = 1.0001 * pressures[0]
    working_gas.SP = working_gas.s, pressures[1]
    working_gas.equilibrate("SP")
    densities[1] = working_gas.density

    # calculate sound speed
    sound_speed = np.sqrt(np.diff(pressures)/np.diff(densities))[0]

    return quant(sound_speed, "m/s")


def calculate_reflected_shock_state(
        initial_temperature,
        initial_pressure,
        species_dict,
        mechanism,
        unit_registry=_U,
        use_multiprocessing=False
):
    """
    Calculates the thermodynamic and chemical state of a reflected shock
    using customized sdtoolbox functions.

    Parameters
    ----------
    initial_temperature : pint.Quantity or Tuple[float, str]
        Mixture initial temperature
    initial_pressure : pint.Quantity or Tuple[float, str]
        Mixture initial pressure
    species_dict : dict
        Dictionary of initial reactant mixture
    mechanism : str
        Mechanism to use for chemical calculations, e.g. ``gri30.cti``
    unit_registry : pint.UnitRegistry, optional
        Pint unit registry
    use_multiprocessing : bool, optional
        True to use multiprocessing for CJ state calculation, which is faster
        but requires the function to be run from ``__main__``

    Returns
    -------
    dict
        Dictionary containing keys

        * ``reflected`` and
        * ``cj``.

        Each of these keys corresponds to a dictionary containing the keys

        * ``speed``, indicating the related wave speed, and
        * ``state``, which is a Cantera gas object at the specified state.
    """
    quant = unit_registry.Quantity

    # define gas objects
    initial_gas = ct.Solution(mechanism)
    reflected_gas = ct.Solution(mechanism)

    # define gas states
    initial_temperature = initial_temperature.to("K").magnitude
    initial_pressure = initial_pressure.to("Pa").magnitude

    initial_gas.TPX = [
        initial_temperature,
        initial_pressure,
        species_dict
    ]
    reflected_gas.TPX = [
        initial_temperature,
        initial_pressure,
        species_dict
    ]

    # get CJ state
    cj_calcs = sd.Detonation.cj_speed(
        initial_pressure,
        initial_temperature,
        species_dict,
        mechanism,
        return_state=True,
        use_multiprocessing=use_multiprocessing
    )

    # get reflected state
    [_,
     reflected_speed,
     reflected_gas] = sd.Reflection.reflect(
        initial_gas,
        cj_calcs["cj state"],
        reflected_gas,
        cj_calcs["cj speed"]
    )

    return {
        "reflected": {
            "speed": quant(
                reflected_speed,
                "m/s"
            ),
            "state": reflected_gas
        },
        "cj": {
            "speed": quant(
                cj_calcs["cj speed"],
                "m/s"),
            "state": cj_calcs["cj state"]
        }
    }


def find_mechanisms(
        return_directory=False
):
    """
    Figure out which mechanisms the local cantera install has access to.

    Parameters
    ----------
    return_directory : bool, optional
        Whether or not to return the location of the mechanism files as well
        as its contents. Defaults to ``False``.

    Returns
    -------
    set or tuple
        Set of available mechanisms in the cantera data directory. If
        `return_directory` is set to True, a tuple is returned where the first
        item is the set of available mechanisms, and the second is the location
        of the cantera data directory.
    """
    mechanism_path = os.path.join(
        os.path.split(os.path.abspath(ct.__file__))[0],
        "data"
    )

    available = {item for item in os.listdir(mechanism_path) if
                 (".cti" in item) or (".xml" in item)}

    if return_directory:
        return available, mechanism_path
    else:
        return available
