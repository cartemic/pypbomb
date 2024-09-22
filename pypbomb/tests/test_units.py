"""
PURPOSE:
    Unit tests for accessories.py

CREATED BY:
    Mick Carter
    Oregon State University
    CIRE and Propulsion Lab
    cartemic@oregonstate.edu
"""

import numpy as np
import pint
import pytest

from .. import units

_U = pint.UnitRegistry()
_Q = _U.Quantity


class TestCheckPintQuantity:
    ureg = pint.UnitRegistry()
    quant = ureg.Quantity

    def test_good_input(self):
        lengths = [6.3, -8]
        for length in lengths:
            assert units.check_pint_quantity(self.quant(length, "inch"), "length")

    def test_bad_dimension_type(self):
        bad_dimension_type = "walrus"
        error_str = bad_dimension_type + " not a supported dimension type"
        with pytest.raises(ValueError, match=error_str):
            units.check_pint_quantity(self.quant(3, "degC"), bad_dimension_type)

    # noinspection PyTypeChecker
    @staticmethod
    def test_non_pint_quantity():
        with pytest.raises(ValueError, match="Non-pint quantity"):
            units.check_pint_quantity(7, "length")

    # noinspection SpellCheckingInspection
    def test_non_numeric_quantity(self):
        with pytest.raises(ValueError, match="Non-numeric pint quantity"):
            units.check_pint_quantity(self.quant("asdf", "inch"), "length")

    def test_ensure_positive_with_negative_magnitude(self):
        with pytest.raises(ValueError, match="Input value < 0"):
            units.check_pint_quantity(self.quant(-4, "in"), "length", ensure_positive=True)

    def test_incorrect_dimensionality(self):
        error_str = (
            self.ureg.degC.dimensionality.__str__().strip("[]")
            + " is not "
            + self.ureg.meter.dimensionality.__str__().strip("[]")
        )
        with pytest.raises(ValueError, match=error_str):
            units.check_pint_quantity(self.quant(19.2, "degC"), "length")


class TestParseQuantInput:
    def test_magnitude(self):
        magnitude = 700
        test_units = "degree_Celsius"
        assert units.parse_quant_input((magnitude, test_units), _U).magnitude == magnitude

    def test_units(self):
        magnitude = 700
        test_units = "degree_Celsius"
        assert units.parse_quant_input((magnitude, test_units), _U).units.format_babel() == test_units

    def test_unit_registry(self):
        magnitude = 700
        test_units = "degree_Celsius"
        new_ureg = pint.UnitRegistry()
        test = units.parse_quant_input((magnitude, test_units), new_ureg)
        assert id(test._REGISTRY) == id(new_ureg)

    def test_pass_quantity(self):
        good_quant = _Q(9, "degC")
        assert units.parse_quant_input(good_quant, _U) == good_quant

    def test_change_registry(self):
        new_ureg = pint.UnitRegistry()
        good_quant = _Q(9, "degC")
        msg = "Cannot operate with Quantity and Quantity of different " "registries."
        with pytest.raises(ValueError, match=msg):
            units.parse_quant_input(good_quant, new_ureg) + good_quant

    def test_bad_input(self):
        bad_inputs = ((0, "psi", "lol"), (0,), "neat")
        checks = np.zeros(len(bad_inputs)).astype(bool)
        for i, bad in enumerate(bad_inputs):
            msg = "Bad quantity input: {0}".format(bad)
            try:
                units.parse_quant_input(bad, _U)
                assert False
            except ValueError as e:
                checks[i] = str(e) == msg
