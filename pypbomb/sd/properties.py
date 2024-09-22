from cantera import Solution


def set_equilibrium(
    working_gas: Solution,
    density: float,
    temperature: float,
):
    """
    Set ``working_gas`` to equilibrium pressure and enthalpy given temperature and density

    :param working_gas: Working gas mixture used for calculations (MUTATED)
    :param density: Mixture density (kg/m^3)
    :param temperature: Mixture temperature in (K)
    :return:
    """
    working_gas.TD = temperature, density
    working_gas.equilibrate("TV")
