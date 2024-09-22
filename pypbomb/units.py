"""
Extra unit management tools for pint quantities
"""

import pint
from numpy import ndarray


def check_pint_quantity(quantity, dimension_type, ensure_positive=False):
    """
    Checks to make sure that a quantity is an instance of a pint quantity, and
    that it has the correct units. Currently supported dimension types:

        * ``"length"``
        * ``"area"``
        * ``"volume"``
        * ``"temperature"``
        * ``"pressure"``
        * ``"velocity"``
        * ``"density"``

    Parameters
    ----------
    quantity : pint.Quantity
        Pint quantity which is to be checked for dimensionality
    dimension_type : str
        Dimensionality that quantity should have
    ensure_positive : bool
        Determines whether the magnitude of the pint quantity will be checked
        for positivity

    Returns
    -------
    bool
        True if no errors are raised
    """
    if hasattr(quantity, "magnitude") and isinstance(quantity.magnitude, ndarray):
        quantity = quantity[0]

    ureg = pint.UnitRegistry()
    units = {
        "length": ureg.meter.dimensionality.__str__(),
        "area": (ureg.meter**2).dimensionality.__str__(),
        "volume": (ureg.meter**3).dimensionality.__str__(),
        "temperature": ureg.degC.dimensionality.__str__(),
        "pressure": ureg.psi.dimensionality.__str__(),
        "velocity": (ureg.meter / ureg.second).dimensionality.__str__(),
        "density": (ureg.kg / ureg.meter**3).dimensionality.__str__(),
    }

    if dimension_type not in units:
        raise ValueError(dimension_type + " not a supported dimension type")

    try:
        actual_dimension_type = quantity.dimensionality.__str__()
    except AttributeError:
        raise ValueError("Non-pint quantity")

    try:
        float(quantity.magnitude)
    except ValueError:
        raise ValueError("Non-numeric pint quantity")

    if ensure_positive:
        if quantity.to_base_units().magnitude < 0:
            raise ValueError("Input value < 0")

    if units[dimension_type] != actual_dimension_type:
        raise ValueError(
            actual_dimension_type.replace("[", "").replace("]", "")
            + " is not "
            + units[dimension_type].replace("[", "").replace("]", "")
        )

    return True


def parse_quant_input(quant_input, unit_registry):
    """
    Converts a tuple of ``(magnitude, "units")`` to a pint quantity or
    converts a pint quantity to the local registry.

    Parameters
    ----------
    quant_input : pint.Quantity or tuple
        Iterable or quantity to be parsed
    unit_registry : pint.UnitRegistry
        Unit registry to be used for pint quantities

    Returns
    -------
    pint.Quantity
        Input as a pint quantity
    """
    if hasattr(quant_input, "magnitude"):
        return unit_registry.Quantity(quant_input.magnitude, quant_input.units.format_babel())
    elif hasattr(quant_input, "__iter__") and len(quant_input) == 2:
        return unit_registry.Quantity(quant_input[0], quant_input[1])
    else:
        raise ValueError("Bad quantity input: {0}".format(quant_input))
