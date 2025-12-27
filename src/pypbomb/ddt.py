"""
Functions relating to the deflagration-to-detonation transition (DDT).
"""

import typing

import cantera as ct
import numpy as np

from pypbomb import _validate, thermochem
from pypbomb._types import SpeciesDefinition


def ideal_blockage_diameter(pipe_id: float, br: float) -> float:
    """
    Calculate the diameter of a Shchelkin spiral corresponding to a given blockage ratio inside a pipe of given inner
    diameter.

    :param pipe_id: Inner diameter of the pipe used for the detonation tube (m)
    :param br: Blockage Ratio: ratio of blocked area to total cross-sectional area, :math:`0 < BR < 1`
    :return: Shchelkin spiral diameter (m)
    """
    _validate.length(pipe_id)
    if not 0 < br < 1:
        raise ValueError(f"Invalid blockage ratio: {br}. Must be on 0 < `blockage_ratio` < 1.")

    return pipe_id / 2 * (1 - np.sqrt(1 - br))


def blockage_ratio(pipe_id: float, blockage_diameter: float) -> float:
    """
    Calculate the blockage ratio of a Shchelkin spiral within a detonation tube.

    :param pipe_id: Inner diameter of the detonation tube (m)
    :param blockage_diameter: Outer diameter of the blockage used to create the Shchelkin spiral (m)
    :return: Ratio of blocked to open area (between 0 and 1)
    """
    _validate.length(pipe_id)
    _validate.length(blockage_diameter)
    if blockage_diameter > pipe_id:
        raise ValueError("Blockage diameter cannot be greater than pipe inner diameter")

    return 1 - (1 - 2 * blockage_diameter / pipe_id) ** 2


def runup_distance(
    pipe_id: float,
    br: float,
    initial_temperature: float,
    initial_pressure: float,
    species: SpeciesDefinition,
    mechanism: str,
    phase_specification: str = "",
) -> float:
    """
    Calculate the runup distance needed for a detonation to develop from a deflagration for a given blockage ratio, tube
    diameter, and mixture. This is accomplished using equations collected by Ciccarelli and Dorofeev [1] for blockage
    ratios <= 0.75. If the desired blockage ratio is less than 0.3, the mixture viscosity is needed, and the
    `phase_specification` parameter may be necessary depending on the mechanism.

    [1] G. Ciccarelli and S. Dorofeev, *Flame acceleration and transition to detonation in ducts*, Progress in Energy
    and Combustion Science, vol. 34, no. 4, pp. 499–550, Aug. 2008. https://doi.org/10.1016/j.pecs.2007.11.002

    :param pipe_id: Inner diameter of the detonation tube (m)
    :param br: Blockage Ratio: ratio of blocked area to total cross-sectional area, :math:`0 < BR <= 0.75`
    :param initial_temperature: Mixture initial temperature (K)
    :param initial_pressure: Mixture initial pressure (Pa)
    :param species: Species definition for cantera
    :param mechanism: String of mechanism to use, e.g. ``gri30.yaml``
    :param phase_specification: Phase specification for cantera solution
    :return: Predicted DDT distance (m)
    """
    _validate.length(pipe_id)
    _validate.temperature(initial_temperature)
    _validate.pressure(initial_pressure)
    if not 0 < br <= 0.75:
        raise ValueError("Blockage ratio outside of correlation range")

    lfs = thermochem.laminar_flame_speed(
        initial_temperature=initial_temperature,
        initial_pressure=initial_pressure,
        species=species,
        mechanism=mechanism,
        phase_specification=phase_specification,
    )

    # calculate density ratio across the deflagration assuming adiabatic flame
    working_gas = ct.Solution(mechanism, phase_specification)
    working_gas.TPX = initial_temperature, initial_pressure, species
    density_reactants = working_gas.density
    working_gas.equilibrate("HP")
    density_products = working_gas.density
    density_ratio = density_reactants / density_products

    # sound speed in products at adiabatic flame temperature
    sound_speed = thermochem.equilibrium_sound_speed(
        temperature=working_gas.T,
        pressure=working_gas.P,
        species=species,
        mechanism=mechanism,
        phase_specification=phase_specification,
    )

    # use appropriate equation to calculate runup distance
    if 0.3 <= br <= 0.75:
        return _runup_eq4_4(
            pipe_id=pipe_id,
            br=br,
            lfs=lfs,
            sound_speed=sound_speed,
            density_ratio=density_ratio,
        )
    elif 0.1 >= br:
        return _runup_eq4_1(
            working_gas=working_gas,
            pipe_id=pipe_id,
            initial_temperature=initial_temperature,
            initial_pressure=initial_pressure,
            species=species,
            br=br,
            lfs=lfs,
            sound_speed=sound_speed,
            density_ratio=density_ratio,
        )
    else:
        high_br_runup = _runup_eq4_1(
            working_gas=working_gas,
            pipe_id=pipe_id,
            initial_temperature=initial_temperature,
            initial_pressure=initial_pressure,
            species=species,
            br=br,
            lfs=lfs,
            sound_speed=sound_speed,
            density_ratio=density_ratio,
        )
        low_br_runup = _runup_eq4_4(
            pipe_id=pipe_id,
            br=br,
            lfs=lfs,
            sound_speed=sound_speed,
            density_ratio=density_ratio,
        )
        interp_distances = np.array([high_br_runup, low_br_runup])
        return typing.cast(float, np.interp(br, np.array([0.1, 0.3]), interp_distances))


def _runup_eq4_1(
    working_gas: ct.Solution,
    pipe_id: float,
    initial_temperature: float,
    initial_pressure: float,
    species: SpeciesDefinition,
    br: float,
    lfs: float,
    sound_speed: float,
    density_ratio: float,
) -> float:
    """
    Calculate runup distance for blockage ratios <= 0.1 using equation 4.1 from G. Ciccarelli and S. Dorofeev, “Flame
    acceleration and transition to detonation in ducts,” Progress in Energy and Combustion Science, vol. 34, no. 4,
    pp. 499–550, Aug. 2008.

    :param working_gas: Working gas mixture used for calculations (MUTATED)
    :param pipe_id: Inner diameter of the detonation tube (m)
    :param initial_temperature: Mixture initial temperature (K)
    :param initial_pressure: Mixture initial pressure (Pa)
    :param species: Species definition for cantera
    :param br: Blockage Ratio: ratio of blocked area to total cross-sectional area, :math:`0 < BR <= 0.1`
    :param lfs: Laminar flame speed (m/s)
    :param sound_speed: Sound speed in products at adiabatic flame temperature (m/s)
    :param density_ratio: Deflagration reactant/product density ratio
    :return: Predicted DDT distance (m)
    """
    # define constants
    kappa = 0.4
    k = 5.5
    c = 0.2
    m = -0.18
    eta = 2.1

    # calculate laminar flame thickness, delta
    working_gas.TPX = initial_temperature, initial_pressure, species
    rho = working_gas.density_mass
    mu = working_gas.viscosity
    nu = mu / rho
    delta = nu / lfs

    # calculate gamma
    gamma = (sound_speed / (eta * (density_ratio - 1) ** 2 * lfs) * (delta / pipe_id) ** (1.0 / 3)) ** (
        1 / (2 * m + 7.0 / 3)
    )

    # calculate runup distance
    d_over_h = 2.0 / (1 - np.sqrt(1 - br))
    return gamma / c * (1 / kappa * np.log(gamma * d_over_h) + k) * pipe_id


def _runup_eq4_4(
    pipe_id: float,
    br: float,
    lfs: float,
    sound_speed: float,
    density_ratio: float,
):
    """
    Calculate runup for blockage ratios between 0.3 and 0.75 using equation 4.4 in G. Ciccarelli and S. Dorofeev, “Flame
    acceleration and transition to detonation in ducts,” Progress in Energy and Combustion Science, vol. 34, no. 4, pp.
    499–550, Aug. 2008.

    :param pipe_id: Inner diameter of the detonation tube (m)
    :param br: Blockage Ratio: ratio of blocked area to total cross-sectional area, :math:`0.3 <= BR <= 0.75`
    :param lfs: Laminar flame speed (m/s)
    :param sound_speed: Sound speed in products at adiabatic flame temperature (m/s)
    :param density_ratio: Deflagration reactant/product density ratio
    :return: Predicted DDT distance (m)
    """
    # define constants
    a = 2.0
    b = 1.5

    # calculate left and right hand sides of eq 4.4
    lhs = 2 * 10 * lfs * (density_ratio - 1) / (sound_speed * pipe_id)
    rhs = a * (1 - br) / (1 + b * br)

    return rhs / lhs
