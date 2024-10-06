"""
Functions for performing thermochemical calculations using ``cantera`` and ``pypbomb.sd``.
"""

import os
from dataclasses import dataclass

import cantera as ct
import numpy as np

from pypbomb import _validate, sd
from pypbomb._types import SpeciesDefinition
from pypbomb.sd.cj import CjShock
from pypbomb.sd.reflect import ReflectedShock


def laminar_flame_speed(
    initial_temperature: float,
    initial_pressure: float,
    species: SpeciesDefinition,
    mechanism: str,
    phase_specification: str = "",
) -> float:
    """
    Calculates the laminar flame speed of a given gas mixture using Cantera.

    :param initial_temperature: Mixture initial temperature (K)
    :param initial_pressure: Mixture initial pressure (Pa)
    :param species: Species definition for cantera
    :param mechanism: String of mechanism to use, e.g. ``gri30.yaml``
    :param phase_specification: Phase specification for cantera solution
    :return: Laminar flame speed (m/s)
    """
    _validate.temperature(initial_temperature)
    _validate.pressure(initial_pressure)

    gas = ct.Solution(mechanism, phase_specification)
    gas.TPX = initial_temperature, initial_pressure, species

    # find laminar flame speed
    flame = ct.FreeFlame(gas)
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
    flame.solve(loglevel=0)

    return flame.velocity[0]


def equilibrium_sound_speed(
    temperature: float,
    pressure: float,
    species: SpeciesDefinition,
    mechanism: str,
    phase_specification: str = "",
):
    """
    Calculate the equilibrium speed of sound in a mixture

    :param temperature: Mixture temperature (K)
    :param pressure: Mixture pressure (Pa)
    :param species: Species definition for cantera
    :param mechanism: Desired chemical mechanism
    :param phase_specification: Phase specification for cantera solution
    :return: Local speed of sound in given mixture at chemical equilibrium
    """
    _validate.temperature(temperature)
    _validate.pressure(pressure)

    gas = ct.Solution(mechanism, phase_specification)
    gas.TPX = [temperature, pressure, species]

    pressures = np.zeros(2)
    densities = np.zeros(2)

    # equilibrate gas at input conditions and collect pressure, density
    gas.equilibrate("TP")
    pressures[0] = gas.P
    densities[0] = gas.density

    # perturb pressure and equilibrate with constant P, s to get dp/drho|s
    pressures[1] = 1.0001 * pressures[0]
    gas.SP = gas.s, pressures[1]
    gas.equilibrate("SP")
    densities[1] = gas.density

    return np.sqrt(np.diff(pressures) / np.diff(densities))[0]


@dataclass(frozen=True)
class ReflectedCjShock:
    cj: CjShock
    reflected: ReflectedShock


def reflected_cj_shock(
    initial_temperature: float,
    initial_pressure: float,
    species: SpeciesDefinition,
    mechanism: str,
    parallelize_cj_calc: bool = False,
) -> ReflectedCjShock:
    """
    Calculates the CJ state along with the state after shock reflection using customized sdtoolbox functions.

    :param initial_temperature: Mixture initial temperature (K)
    :param initial_pressure: Mixture initial pressure (Pa)
    :param species: Species definition for cantera
    :param mechanism: Mechanism to use for chemical calculations, e.g. ``gri30.yaml``
    :param parallelize_cj_calc: Use multiprocessing for CJ state calculation, which is faster but requires the function
        to be run from ``__main__``, which may not behave well when used via Jupyter
    :return:
    """
    _validate.temperature(initial_temperature)
    _validate.pressure(initial_pressure)

    initial_gas = ct.Solution(mechanism)
    initial_gas.TPX = [initial_temperature, initial_pressure, species]

    cj = sd.cj.shock(
        initial_pressure=initial_pressure,
        initial_temperature=initial_temperature,
        species=species,
        mechanism=mechanism,
        parallelize=parallelize_cj_calc,
    )
    reflected = sd.reflect.shock(
        mechanism=mechanism,
        initial_state_gas=initial_gas,
        post_shock_gas=cj.state,
        incident_shock_speed=cj.speed,
    )

    return ReflectedCjShock(cj=cj, reflected=reflected)


@dataclass(frozen=True)
class FoundMechanisms:
    mechanisms: set[str]
    path: str


def find_mechanisms(filetypes: tuple[str] = (".cti", ".xml", ".yaml")) -> FoundMechanisms:
    """
    Figure out which mechanisms the local cantera install has access to.

    :return:
    """
    path = os.path.join(os.path.split(os.path.abspath(ct.__file__))[0], "data")
    mechanisms = {item for item in os.listdir(path) if item.endswith(filetypes)}

    return FoundMechanisms(mechanisms=mechanisms, path=path)
