# -*- coding: utf-8 -*-

import os
from math import sqrt

import cantera as ct
import numpy as np
import pandas as pd
import pytest
from mock import patch

from .. import thermochem, tube


def compare(manual, tested):
    for key, value in manual.items():
        try:
            test_value = tested[key].to(value.units.format_babel()).magnitude
            value = value.magnitude
        except AttributeError:
            test_value = tested[key]

        assert np.all(np.abs(test_value - value) / value < 1e-4)


# noinspection SpellCheckingInspection
class TestBolt:
    thread_size = "1/4-28"
    thread_class = "2"
    plate_max_tensile = tube._Q(30, "ksi")
    engagement_length = tube._Q(0.5, "in")
    bolt_max_tensile = tube._Q(80, "ksi")

    def test_calculate_sf(self):
        max_pressure = tube._Q(1631.7, "psi")
        window_area = tube._Q(5.75 * 2.5, "in^2")
        num_bolts = 20
        thread_size = "1/4-28"
        thread_class = "2"
        bolt_max_tensile = tube._Q(120, "ksi")
        plate_max_tensile = (30, "ksi")  # parse check
        engagement_length = tube._Q(0.5, "in")

        hand_calc = {
            "bolt": 3.606968028,
            "plate": 7.969517321,
        }

        test_values = tube.Bolt.calculate_safety_factors(
            max_pressure,
            window_area,
            num_bolts,
            thread_size,
            thread_class,
            bolt_max_tensile,
            plate_max_tensile,
            engagement_length,
            tube._U,
        )

        compare(hand_calc, test_values)

    def test_calculate_sf_vectorized(self):
        max_pressure = tube._Q(1631.7, "psi")
        window_area = tube._Q(5.75 * 2.5, "in^2")
        num_bolts = 20 * np.ones(12)
        thread_size = "1/4-28"
        thread_class = "2"
        bolt_max_tensile = tube._Q(120, "ksi")
        plate_max_tensile = tube._Q(30, "ksi")
        engagement_length = tube._Q(0.5, "in")

        hand_calc = {
            "bolt": 3.606968028,
            "plate": 7.969517321,
        }

        test_values = tube.Bolt.calculate_safety_factors(
            max_pressure,
            window_area,
            num_bolts,
            thread_size,
            thread_class,
            bolt_max_tensile,
            plate_max_tensile,
            engagement_length,
            tube._U,
        )

        compare(hand_calc, test_values)

    def test_calculate_stress_areas_over_100ksi(self):
        # test bolt > 100ksi
        bolt_max_tensile = tube._Q(120, "ksi")
        hand_calc = {
            "screw area": tube._Q(0.034934049, "in^2"),
            "plate area": tube._Q(0.308744082, "in^2"),
            "minimum engagement": tube._Q(0.452595544, "in"),
        }
        test_areas = tube.Bolt.calculate_stress_areas(
            self.thread_size,
            self.thread_class,
            bolt_max_tensile,
            self.plate_max_tensile,
            self.engagement_length,
            tube._U,
        )
        compare(hand_calc, test_areas)

    def test_calculate_stress_areas_under_100ksi(self):
        # test bolt < 100ksi
        bolt_max_tensile = tube._Q(80, "ksi")
        hand_calc = {
            "screw area": tube._Q(0.036374073, "in^2"),
            "minimum engagement": tube._Q(0.314168053, "in"),
        }
        test_areas = tube.Bolt.calculate_stress_areas(
            self.thread_size,
            self.thread_class,
            bolt_max_tensile,
            self.plate_max_tensile,
            self.engagement_length,
            tube._U,
        )
        compare(hand_calc, test_areas)

    def test_calculate_stress_areas_length_too_short(self):
        # ensure warning when engagement length < minimum
        engagement_length = tube._Q(0.05, "in")
        with pytest.warns(
            Warning,
            match="Screws fail in shear, not tension."
            + " Plate may be damaged."
            + " Consider increasing bolt engagement length",
        ):
            tube.Bolt.calculate_stress_areas(
                self.thread_size,
                self.thread_class,
                self.bolt_max_tensile,
                self.plate_max_tensile,
                engagement_length,
                tube._U,
            )

    @staticmethod
    def test_get_thread_property():
        # good input
        good_args = [
            ["pitch diameter max", "1/4-20", "2B", tube._U],
            ["pitch diameter max", "1/4-20", "2A", tube._U],
        ]
        good_results = [0.2248, 0.2164]
        for args, result in zip(good_args, good_results):
            test_result = tube.Bolt.get_thread_property(*args)

            assert abs(test_result.magnitude - result) < 1e-7

    @staticmethod
    def test_get_thread_property_not_in_dataframe():
        # property not in dataframe
        bad_property = "jello"
        bad_message = (
            "Thread property '"
            + bad_property
            + "' not found. Available specs: "
            + "'"
            + "', '".join(tube.THREAD_SPECS["internal"].keys())
            + "'"
        )
        with pytest.raises(KeyError, match=bad_message):
            tube.Bolt.get_thread_property(bad_property, "1/4-20", "2B", tube._U)

    @staticmethod
    def test_get_thread_property_invalid_class():
        bad_args = [
            ["pitch diameter max", "1/4-20", "2F", tube._U],
            ["pitch diameter max", "1/4-20", "6A", tube._U],
        ]
        for args in bad_args:
            with pytest.raises(ValueError, match="bad thread class"):
                tube.Bolt.get_thread_property(*args)

    @staticmethod
    def test_get_thread_property_invalid_size():
        bad_size = "1290-33"
        bad_message = "Thread size '" + bad_size + "' not found"
        with pytest.raises(KeyError, match=bad_message):
            tube.Bolt.get_thread_property("pitch diameter max", bad_size, "2B", tube._U)


# noinspection SpellCheckingInspection
class TestDDT:
    diameter = tube._Q(5.76, tube._U.inch)

    # use a unit diameter to match diameter-specific values from plot
    tube_diameter = tube._Q(1, "meter")

    # define gas mixture and relevant pint quantities
    mechanism = "gri30.cti"
    gas = ct.Solution(mechanism)
    gas.TP = 300, 101325
    initial_temperature = tube._Q(gas.T, "K")
    initial_pressure = tube._Q(gas.P, "Pa")
    gas.set_equivalence_ratio(1, "CH4", {"O2": 1, "N2": 3.76})
    species_dict = gas.mole_fraction_dict()

    def test_calculate_spiral_diameter(self):
        # define a blockage ratio
        test_blockage_ratio = 0.44

        # define expected result and actual result
        expected_spiral_diameter = self.diameter / 2 * (1 - sqrt(1 - test_blockage_ratio))
        result = tube.DDT.calculate_blockage_diameter(self.diameter, test_blockage_ratio)

        # ensure good output
        assert expected_spiral_diameter == result.to(tube._U.inch)

    def test_calculate_spiral_diameter_non_numeric_br(self):
        # ensure proper handling with non-numeric blockage ratio
        with pytest.raises(ValueError, match="Non-numeric blockage ratio."):
            # noinspection PyTypeChecker
            tube.DDT.calculate_blockage_diameter(self.diameter, "doompity doo")

    def test_calculate_spiral_diameter_bad_br(self):
        # ensure proper handling with blockage ratio outside allowable limits
        bad_blockage_ratios = [-35.124, 0, 1, 120.34]
        for ratio in bad_blockage_ratios:
            with pytest.raises(ValueError, match="Blockage ratio outside of 0<BR<1"):
                tube.DDT.calculate_blockage_diameter(self.diameter, ratio)

    @staticmethod
    def test_calculate_blockage_ratio():
        hand_calc_blockage_ratio = 0.444242326717679

        # check for expected result with good input
        test_result = tube.DDT.calculate_blockage_ratio(tube._Q(3.438, tube._U.inch), tube._Q(7.0 / 16.0, tube._U.inch))
        assert (test_result - hand_calc_blockage_ratio) < 1e-8

    @staticmethod
    def test_calculate_blockage_ratio_0_tube_diameter():
        with pytest.raises(ValueError, match="tube ID cannot be 0"):
            tube.DDT.calculate_blockage_ratio(tube._Q(0, tube._U.inch), tube._Q(7.0 / 16.0, tube._U.inch))

    def test_calculate_blockage_ratio_blockage_gt_tube(self):
        # check for correct handling when blockage diameter >= tube diameter
        with pytest.raises(ValueError, match="blockage diameter >= tube diameter"):
            tube.DDT.calculate_blockage_ratio(tube._Q(1, tube._U.inch), tube._Q(3, tube._U.inch))

    def test_calculate_ddt_run_up_bad_blockage_ratio(self):
        bad_blockages = [-4.0, 0, 1]
        for blockage_ratio in bad_blockages:
            with pytest.raises(ValueError, match="Blockage ratio outside of correlation range"):
                tube.DDT.calculate_run_up(
                    blockage_ratio,
                    self.tube_diameter,
                    self.initial_temperature,
                    self.initial_pressure,
                    self.species_dict,
                    self.mechanism,
                    tube._U,
                )

    def test_calculate_ddt_run_up(self):
        # define good blockage ratios and expected result from each
        good_blockages = [0.1, 0.2, 0.3, 0.75]
        good_results = [
            48.51385390428211,
            29.24433249370277,
            18.136020151133494,
            4.76070528967254,
        ]

        # test with good inputs
        for blockage_ratio, result in zip(good_blockages, good_results):
            test_runup = tube.DDT.calculate_run_up(
                blockage_ratio,
                self.tube_diameter,
                self.initial_temperature,
                self.initial_pressure,
                self.species_dict,
                self.mechanism,
                tube._U,
                phase_specification="gri30_mix",
            )

            assert 0.5 * result <= test_runup.magnitude <= 1.5 * result


class TestWindow:
    width = tube._Q(50, tube._U.mm).to(tube._U.inch)
    length = tube._Q(20, tube._U.mm).to(tube._U.mile)
    pressure = tube._Q(1, tube._U.atm).to(tube._U.torr)
    thickness = tube._Q(1.2, tube._U.mm).to(tube._U.furlong)
    rupture_modulus = tube._Q(5300, tube._U.psi).to(tube._U.mmHg)

    def test_safety_factor(self):
        desired_safety_factor = 4

        test_sf = tube.Window.safety_factor(
            self.length, self.width, self.thickness, self.pressure, self.rupture_modulus
        )

        assert abs(test_sf - desired_safety_factor) / test_sf < 0.01

    def test_safety_factor_vectorized(self):
        desired_safety_factor = 4

        test_sf = tube.Window.safety_factor(
            self.length,
            (self.width.magnitude, self.width.units.format_babel()),  # parse
            self.thickness,
            self.pressure,
            self.rupture_modulus * np.ones(7),  # vectorization check
        )

        assert np.all(np.abs(test_sf - desired_safety_factor) / test_sf < 0.01)

    def test_minimum_thickness_sf_less_than_1(self):
        # safety factor < 1
        safety_factor = [0.25, -7]
        for factor in safety_factor:
            with pytest.raises(ValueError, match="Window safety factor < 1"):
                tube.Window.minimum_thickness(
                    self.length,
                    self.width,
                    factor,
                    self.pressure,
                    self.rupture_modulus,
                    tube._U,
                )

    def test_minimum_thickness_sf_non_numeric(self):
        safety_factor = "BruceCampbell"
        with pytest.raises(TypeError, match="Non-numeric window safety factor"):
            # noinspection PyTypeChecker
            tube.Window.minimum_thickness(
                self.length,
                self.width,
                safety_factor,
                self.pressure,
                self.rupture_modulus,
                tube._U,
            )

    def test_minimum_thickness(self):
        safety_factor = 4
        test_thickness = tube.Window.minimum_thickness(
            self.length,
            self.width,
            safety_factor,
            self.pressure,
            self.rupture_modulus,
            tube._U,
        )
        test_thickness = test_thickness.to(self.thickness.units).magnitude
        desired_thickness = self.thickness.magnitude
        assert abs(test_thickness - desired_thickness) / test_thickness < 0.01

    def test_minimum_thickness_vectorized(self):
        safety_factor = 4
        test_thickness = tube.Window.minimum_thickness(
            self.length * np.ones(9),  # vectorization check
            (self.width.magnitude, self.width.units.format_babel()),  # parse
            safety_factor,
            self.pressure,
            self.rupture_modulus,
            tube._U,
        )
        test_thickness = test_thickness.to(self.thickness.units).magnitude
        desired_thickness = self.thickness.magnitude
        assert np.all(np.abs(test_thickness - desired_thickness) / test_thickness < 0.01)


class TestTube:
    def test__get_pipe_stress_limits(self):
        unit_registry = tube._U
        test_limits_weld = tube.Tube._get_pipe_stress_limits("316", True, unit_registry)
        correct_316_weld = [
            16,
            16,
            13.8,
            12.4,
            11.4,
            10.6,
            10.1,
            9.8,
            9.6,
            9.5,
            9.4,
            9.3,
            9.2,
            9.1,
            9,
            8.9,
            8.8,
            7.9,
            6.3,
        ]
        test_limits_smls = tube.Tube._get_pipe_stress_limits("316", False, unit_registry)
        correct_316_smls = [
            18.8,
            18.8,
            16.2,
            14.6,
            13.4,
            12.5,
            11.8,
            11.6,
            11.3,
            11.2,
            11,
            10.9,
            10.8,
            10.7,
            10.6,
            10.5,
            10.3,
            9.3,
            7.4,
        ]

        test_limits_weld = [t.to("ksi").magnitude for t in test_limits_weld]
        test_limits_smls = [t.to("ksi").magnitude for t in test_limits_smls]

        assert np.allclose(
            np.concatenate((test_limits_weld, test_limits_smls)),
            np.concatenate((correct_316_weld, correct_316_smls)),
        )

    def test_dynamic_load_factor(self):
        quant = tube._Q

        # hand calcs use:
        #     material = "316L"
        #     schedule = "80"
        #     nominal_size = "6"
        tube_id = quant(5.761, "in")
        tube_od = quant(6.625, "in")
        elastic_modulus = quant(200, "GPa")
        density = quant(7.9, "g/cm**3")
        poisson_ratio = 0.28

        # from hand calcs, critical velocity is 1457.44 m/s, giving upper and
        # lower bounds of 1603.188 and 1311.700 m/s
        cj_velocities = [
            quant(1200, "m/s"),  # DLF 1
            quant(1311, "m/s"),  # DLF 1
            quant(1312, "m/s"),  # DLF 4
            quant(1400, "m/s"),  # DLF 4
            quant(1603, "m/s"),  # DLF 4
            quant(1604, "m/s"),  # DLF 2
            quant(2000, "m/s"),  # DLF 2
        ]
        expected_dlf = [1, 1, 4, 4, 4, 2, 2]
        test_dlf = []
        for cj_velocity in cj_velocities:
            test_dlf.append(
                tube.Tube.dynamic_load_factor(
                    tube_id,
                    tube_od,
                    cj_velocity,
                    elastic_modulus,
                    density,
                    poisson_ratio,
                    plus_or_minus=0.1,
                )
            )

        assert all(dlf == good_dlf for dlf, good_dlf in zip(test_dlf, expected_dlf))

    def test_dynamic_load_factor_bad_plus_minus_value(self):
        quant = tube._Q

        # hand calcs use:
        #     material = "316L"
        #     schedule = "80"
        #     nominal_size = "6"
        tube_id = quant(5.761, "in")
        tube_od = quant(6.625, "in")
        elastic_modulus = quant(200, "GPa")
        density = quant(7.9, "g/cm**3")
        poisson_ratio = 0.28
        cj_velocity = quant(2, "m/s")  # unimportant

        # plus_or_minus outside of (0, 1)
        bad_plus_minus = [-1, 0, 1, 2]
        for plus_minus in bad_plus_minus:
            with pytest.raises(ValueError, match="plus_or_minus factor not between 0 and 1"):
                tube.Tube.dynamic_load_factor(
                    tube_id,
                    tube_od,
                    cj_velocity,
                    elastic_modulus,
                    density,
                    poisson_ratio,
                    plus_minus,
                )

    def test_calculate_initial_pressure_without_mp(self):
        unit_registry = tube._U
        quant = tube._Q
        #     material = "316L"
        #     schedule = "80"
        #     nominal_size = "6"
        initial_temperature = quant(300, "K")
        tube_id = quant(5.761, "in")
        tube_od = quant(6.625, "in")
        elastic_modulus = quant(200, "GPa")
        density = quant(7.9, "g/cm**3")
        poisson_ratio = 0.28
        reactant_mixture = dict(H2=1 / 3, O2=2 / 3)
        mechanism = "gri30.cti"

        # the initial pressure should cause the reflected detonation pressure
        # to be equal to the tube's max pressure, accounting for dynamic load
        # factor
        correct_max = 15101889.83212825  # Pa
        calc_initial = tube.Tube.calculate_max_initial_pressure(
            tube_id,
            tube_od,
            initial_temperature,
            reactant_mixture,
            mechanism,
            quant(correct_max, "Pa"),
            elastic_modulus,
            density,
            poisson_ratio,
            unit_registry=unit_registry,
            use_multiprocessing=False,
        )
        test_state = thermochem.calculate_reflected_shock_state(
            initial_temperature,
            calc_initial,
            reactant_mixture,
            mechanism,
            unit_registry,
        )
        dlf = 2.0

        error = abs(correct_max - test_state["reflected"]["state"].P * dlf) / correct_max

        assert error <= 0.0005

    def test_calculate_initial_pressure_with_mp(self):
        unit_registry = tube._U
        quant = tube._Q
        #     material = "316L"
        #     schedule = "80"
        #     nominal_size = "6"
        tube_id = quant(5.761, "in")
        tube_od = quant(6.625, "in")
        initial_temperature = quant(300, "K")
        reactant_mixture = dict(H2=2 / 3, O2=1 / 3)
        mechanism = "gri30.cti"
        elastic_modulus = quant(200, "GPa")
        density = quant(7.9, "g/cm**3")
        poisson_ratio = 0.28

        # the initial pressure should cause the reflected detonation pressure
        # to be equal to the tube's max pressure, accounting for dynamic load
        # factor
        correct_max = 15101889.83212825  # Pa
        calc_initial = tube.Tube.calculate_max_initial_pressure(
            tube_id,
            tube_od,
            initial_temperature,
            reactant_mixture,
            mechanism,
            quant(correct_max, "Pa"),
            elastic_modulus,
            density,
            poisson_ratio,
            unit_registry=unit_registry,
            use_multiprocessing=True,
        )
        test_state = thermochem.calculate_reflected_shock_state(
            initial_temperature,
            calc_initial,
            reactant_mixture,
            mechanism,
            unit_registry,
        )
        dlf = 2.0

        error = abs(correct_max - test_state["reflected"]["state"].P * dlf) / correct_max

        assert error <= 0.0005

    def test_get_available_pipe_schedules(self):
        good = {
            "8": [
                "5s",
                "5",
                "10s",
                "10",
                "20",
                "30",
                "40s",
                "40",
                "60",
                "80s",
                "80",
                "100",
                "120",
                "140",
                "160",
                "XXH",
            ],
            "26": ["10", "20", "40s", "80s"],
        }
        checks = np.zeros(len(good)).astype(bool)
        for i, (size, correct) in enumerate(good.items()):
            checks[i] = tube.Tube.get_available_pipe_schedules(size) == correct

        assert np.all(checks)

    def test_get_available_pipe_schedules_bad_size(self):
        msg = "Invalid pipe size: 2.222. See Tube.available_pipe_sizes."
        with pytest.raises(ValueError, match=msg):
            # noinspection PyTypeChecker
            tube.Tube.get_available_pipe_schedules(2.222)

    def test_get_dimensions(self):
        # use NPS 9
        size = "9"
        good = {
            "40s": dict(
                inner_diameter=tube._Q(8.941, "in"),
                outer_diameter=tube._Q(9.625, "in"),
                wall_thickness=tube._Q(0.342, "in"),
            ),
            "80s": dict(
                inner_diameter=tube._Q(8.625, "in"),
                outer_diameter=tube._Q(9.625, "in"),
                wall_thickness=tube._Q(0.5, "in"),
            ),
        }
        checks = np.zeros(len(good)).astype(bool)
        for i, (schedule, dims) in enumerate(good.items()):
            checks[i] = tube.Tube.get_dimensions(size, schedule) == dims

        assert np.all(checks)

    def test_get_dimensions_bad_schedule(self):
        schedule = "20"
        sizes = ("1/8", "4 1/2")
        checks = np.zeros(len(sizes)).astype(bool)
        for i, size in enumerate(sizes):
            msg = (
                "Schedule {:s} invalid for pipe size {:s}. See "
                "Tube.available_pipe_sizes and "
                "Tube.get_available_pipe_schedules.".format(schedule, size)
            )
            try:
                tube.Tube.get_dimensions(size, schedule)
                assert False
            except ValueError as e:
                checks[i] = str(e) == msg

        assert np.all(checks)

    def test_calculate_max_stress_seamless(self):
        material = "316L"
        initial_temperatures = [
            tube._Q(100, "degF"),
            tube._Q(200, "degF"),
            tube._Q(150, "degF"),
        ]

        # in ksi
        good_stresses = [15.7, 13.3, 14.5]
        for temperature, good_stress in zip(initial_temperatures, good_stresses):
            test_stress = tube.Tube.calculate_max_stress(temperature, material, False, tube._U).to("ksi").magnitude
            assert np.allclose(test_stress, good_stress)

    def test_calculate_max_stress_welded(self):
        material = "304"
        initial_temperatures = [
            tube._Q(650, "degF"),
            tube._Q(700, "degF"),
            tube._Q(675, "degF"),
        ]

        # in ksi
        good_stresses = [9.5, 9.4, 9.45]
        for temperature, good_stress in zip(initial_temperatures, good_stresses):
            test_stress = tube.Tube.calculate_max_stress(temperature, material, True, tube._U).to("ksi").magnitude
            assert np.allclose(test_stress, good_stress)

    def test_calculate_max_stress_non_monotonic(self):
        material = "316"
        temperature = tube._Q(100, "degF")

        def make_non_monotonic(*_, **__):
            df_fake_stress = pd.DataFrame(
                columns=("temperature", "stress"),
                data=np.array([[0, 1, -1], [0, 0, 0]]).T,
            )
            df_fake_stress["temperature"] = df_fake_stress["temperature"].apply(lambda x: tube._Q(x, "degF"))
            df_fake_stress["stress"] = df_fake_stress["stress"].apply(lambda x: tube._Q(x, "ksi"))
            return df_fake_stress.set_index("temperature")["stress"]

        with patch("pypbomb.tube.Tube._get_pipe_stress_limits", new=make_non_monotonic):
            with pytest.raises(
                ValueError,
                match="Stress limits require temperatures to be " + "monotonically increasing",
            ):
                tube.Tube.calculate_max_stress(temperature, material, True, tube._U)

    def test_get_elastic_modulus(self):
        ureg = thermochem._U
        test = tube.Tube.get_elastic_modulus("316L", ureg)
        assert test == ureg.Quantity(200, "GPa")

    def test_get_density(self):
        ureg = thermochem._U
        test = tube.Tube.get_density("316L", ureg)
        assert test == ureg.Quantity(7.9, "g/cm^3")

    def test_get_poisson(self):
        test = tube.Tube.get_poisson("316L")
        assert test == 0.28

    def test__check_material_good(self):
        test = tube.Tube._check_material("304L")
        assert test is None

    def test__check_material_bad(self):
        bad_material = "unobtainium"
        msg = "Material %s not found" % bad_material
        with pytest.raises(ValueError, match=msg):
            tube.Tube._check_material(bad_material)


class TestCollectTubeMaterials:
    def test_good_output(self):
        columns = ["Grade", "Group", "ElasticModulus", "Density", "Poisson"]
        data = [
            ["304", 2.1, 200, 7.8, 0.28],
            ["304H", 2.1, 200, 7.8, 0.28],
            ["316", 2.2, 200, 7.9, 0.28],
            ["316H", 2.2, 200, 7.9, 0.28],
            ["317", 2.2, 200, 7.9, 0.28],
            ["304L", 2.3, 200, 7.8, 0.28],
            ["316L", 2.3, 200, 7.9, 0.28],
        ]
        good_dataframe = pd.DataFrame(data=data, index=None, columns=columns)
        good_dataframe["ElasticModulus"] = [tube._Q(item, "GPa") for item in good_dataframe["ElasticModulus"]]
        good_dataframe["Density"] = [tube._Q(item, "g/cm^3") for item in good_dataframe["Density"]]
        df_test = tube._collect_tube_materials()

        assert df_test.equals(good_dataframe)

    def test_nonexistent_file(self):
        def fake_exists(*_):
            return False

        with patch("os.path.exists", new=fake_exists):
            with pytest.raises(ValueError, match="materials_list.csv does not exist"):
                tube._collect_tube_materials()

    def test_empty_file(self):
        # noinspection PyUnresolvedReferences
        def fake_csv_import(*_):
            raise pd.errors.EmptyDataError

        with patch("pandas.read_csv", new=fake_csv_import):
            # check for proper error handling when file is blank
            with pytest.raises(ValueError, match="materials_list.csv is empty"):
                tube._collect_tube_materials()

    def test_calculate_max_pressure(self):
        safety_factor = 4

        max_stress = tube._Q(18.8, "ksi")
        wall_thickness = tube._Q(0.432, "in")
        tube_od = tube._Q(6.625, "in")
        tube_id = tube_od - 2 * wall_thickness
        mean_diameter = (tube_id + tube_od) / 2

        good_max_pressure = (max_stress * (2 * wall_thickness) / (mean_diameter * safety_factor)).to("Pa").magnitude

        test_max_pressure = (
            tube.Tube.calculate_max_pressure(tube_id, tube_od, max_stress, safety_factor).to("Pa").magnitude
        )

        assert np.allclose(test_max_pressure, good_max_pressure)


class TestCheckMaterial:
    def test_good_material(self):
        material = "316"
        assert tube._check_material(material)

    def test_bad_material(self):
        material = "unobtainium"
        match_string = "%s not found in materials." " Try Flange.available_materials" % material
        try:
            tube._check_material(material)
            assert False
        except ValueError as err:
            err = str(err)
            assert match_string == err


class TestCheckClass:
    def test__check_class_good(self):
        flange_class = "900"
        group = "2.2"
        assert tube.Flange._check_flange_class(flange_class, group)

    def test__check_class_bad(self):
        flange_class = "9001"
        group = "2.2"
        match_string = "{:s} not found in group {:s}.".format(
            str(flange_class),
            str(group),
        )

        try:
            tube.Flange._check_flange_class(flange_class, group)
            assert False
        except ValueError as err:
            err = str(err)
            assert match_string == err


class TestGetFlangeLimitsFromCSV:
    def test_good_output(self):
        flange_class = "900"
        temperature = tube._Q(400, "degC")

        # good pressures manually entered based on values from .csv files
        good_pressures = {"2.1": 85.3, "2.2": 88.3, "2.3": 72.9}

        result = tube._get_flange_limits_from_csv()

        tests = []
        for group, pressure in good_pressures.items():
            # read the imported pressure and compare it to the correct pressure
            test_pressure = result[group].loc[result[group]["Temperature"] == temperature][flange_class].values[0]
            good_pressure = tube._Q(pressure, "bar")
            tests.append(test_pressure == good_pressure)

        assert all(tests)

    def test_bad_material(self):
        file_name = "ASME_B16_5_flange_ratings_group_2_1.csv"
        file_directory = os.path.join(os.path.dirname(os.path.relpath(tube.__file__)), "lookup_data")
        file_location = os.path.join(file_directory, file_name)
        match_string = "\n" + file_location + "not found"

        # bad group is identified by checking the existence of a .csv fil
        def fake_os_path_exists(*_):
            return False

        with patch("os.path.exists", new=fake_os_path_exists):
            try:
                tube._get_flange_limits_from_csv()
            except FileNotFoundError as err:
                err = str(err)
                assert match_string == err

    def test_bad_pressure(self):
        match_string = "\nPressure less than zero."

        # an error should be thrown when a non-temperature column has a
        # negative value
        def fake_read_csv(*_):
            return pd.DataFrame(columns=["Temperature", "900"], data=[[10, -10]])

        with patch("pandas.read_csv", new=fake_read_csv):
            with pytest.raises(ValueError, match=match_string):
                tube._get_flange_limits_from_csv()

    def test_zeroed_non_numeric(self):
        # create test dataframe and a dataframe of what the output should
        # look like after the test dataframe is imported
        test_temperatures = [9, "s", "d"]
        good_temperatures = [9.0, 0.0, 0.0]
        test_pressures = ["a", 3, "f"]
        good_pressures = [0.0, 3.0, 0.0]
        test_dataframe = pd.DataFrame({"Temperature": test_temperatures, "Pressure": test_pressures})
        good_dataframe = pd.DataFrame({"Temperature": good_temperatures, "Pressure": good_pressures})

        # convince the tube to load the test dataframe without having to
        # write a .csv to disk
        def fake_read_csv(*_):
            return test_dataframe

        with patch("pandas.read_csv", new=fake_read_csv):
            tube._get_flange_limits_from_csv()

            tests = []
            for group in set(tube.TUBE_MATERIALS.Group.values):
                # pull the dataframe of flange limits for each material group
                # and use applymap to remove units from all pint quantities.
                # This is done to circumvent an error caused by using
                # np.allclose for an array of pint quantities. Comparison to
                # good_dataframe.values is done to avoid having to care about
                # how the values are distributed by pandas.
                df_test = tube._get_flange_limits_from_csv()[str(group)].applymap(lambda x: x.magnitude)
                test = np.allclose(df_test, good_dataframe.values)
                tests.append(test)

        assert all(tests)


class TestFlange:
    def test_get_class_good_output(self):
        flange_class = "900"
        temperature = tube._Q(400, "degC")

        # good pressures manually entered based on values from .csv files
        good_pressures = {
            "304": (85.3, "bar"),  # 2.1
            "316": tube._Q(88.3, "bar"),  # 2.2
            "304L": tube._Q(72.9, "bar"),  # 2.3
        }

        tests = []
        for material, pressure in good_pressures.items():
            test_class = tube.Flange.get_class(pressure, temperature, material)

            tests.append(test_class == flange_class)

        assert all(tests)

    def test_get_class_bad_temperature(self):
        # bad temperatures for 316L
        test_temperatures = [(-100, "degC"), tube._Q(500, "degC")]
        material = "316L"
        pressure = tube._Q(125, "bar")

        for temperature in test_temperatures:
            with pytest.raises(ValueError, match="Temperature out of range."):
                tube.Flange.get_class(pressure, temperature, material)

    def test_get_class_bad_pressure(self):
        temperature = (350, "degC")
        pressure = (350, "bar")
        material = "316L"
        with pytest.raises(ValueError, match="Pressure out of range."):
            tube.Flange.get_class(pressure, temperature, material)

    def test_get_max_pressure_good(self):
        temps = [tube._Q(200, "degC"), (225, "degC")]
        test_class = "150"
        material = "304"
        good = [13.2, 12.65]  # bar

        results = np.zeros(2)
        for i, t in enumerate(temps):
            p_test = tube.Flange.get_max_pressure(test_class, t, material)
            results[i] = p_test.to("bar").magnitude

        assert np.allclose(results, good)

    def test_get_max_pressure_bad_temperature(self):
        with pytest.raises(ValueError, match="Temperature out of range"):
            tube.Flange.get_max_pressure("150", tube._Q(10000, "degC"), "316")

    def test_get_max_temperature_good(self):
        pressures = [tube._Q(14.8, "bar"), (14.25, "bar")]
        test_class = "150"
        material = "316"
        good = [150, 175]  # degC

        results = np.zeros(2)
        for i, p in enumerate(pressures):
            p_test = tube.Flange.get_max_temperature(test_class, p, material)
            results[i] = p_test.to("degC").magnitude

        assert np.allclose(results, good)

    def test_get_max_temperature_bad_pressure(self):
        with pytest.raises(ValueError, match="Pressure out of range"):
            tube.Flange.get_max_temperature("150", tube._Q(19.1, "bar"), "316")


class TestCheckMaterialList:
    fake_material_groups = {"thing0": "group0", "thing1": "group1"}
    fake_materials = pd.DataFrame(columns=["Group", "Grade"])

    class FakeOpen:
        """
        fake open()
        """

        def __init__(self, *_):
            """
            dummy init statement
            """

        def __enter__(self, *_):
            """
            enter statement that returns a FakeFile
            """
            return self.FakeFile()

        def __exit__(self, *_):
            """
            dummy exit statement
            """
            return None

        # noinspection SpellCheckingInspection,SpellCheckingInspection
        class FakeFile:
            """
            fake file used for FakeOpen
            """

            @staticmethod
            def readline(*_):
                """
                fake file for use with FakeOpen()
                """
                return "ASDF,thing0,thing1\n"

    def test_good_result(self):
        assert tube._check_materials_list()

    def test_no_files_found(self):
        def fake_listdir(*_):
            # returns an empty list to test for files-not-found condition
            return []

        with pytest.raises(FileNotFoundError, match="No files containing 'flange' or 'stress' found"):
            with patch("os.listdir", new=fake_listdir):
                tube._check_materials_list()

    # noinspection SpellCheckingInspection
    def test_no_welded_or_seamless(self):
        # noinspection SpellCheckingInspection
        def fake_listdir(*_):
            # causes welded/seamless warning
            return ["asdfflangegroup0sfh", "asdfflangegroup1asd", "asdfstresswdasdg"]

        with patch("builtins.open", new=self.FakeOpen):
            with patch("os.listdir", new=fake_listdir):
                msg = "asdfstresswdasdg does not indicate whether it is " "welded or seamless"
                with pytest.warns(UserWarning, match=msg):
                    try:
                        tube._check_materials_list()
                    except ValueError:
                        # ignore value error from non-existent data
                        pass

    # noinspection SpellCheckingInspection
    def test_check_materials_list_missing_material(self):
        file_directory = os.path.join("pypbomb", "lookup_data")

        # noinspection SpellCheckingInspection
        def fake_listdir(*_):
            """
            listdir function that should work 100%
            """
            return [
                "asdfflangegroup0sfh",
                "asdfflangegroup1asd",
                "asdfstressweldedasdg",
            ]

        # test_tube._material_groups = self.fake_material_groups
        old_materials = tube.TUBE_MATERIALS.copy()
        tube.TUBE_MATERIALS = pd.DataFrame(
            columns=["Grade", "Group"],
            data=list(
                map(
                    list,
                    zip(
                        self.fake_material_groups.keys(),
                        self.fake_material_groups.values(),
                    ),
                )
            ),
        )

        with patch("builtins.open", new=self.FakeOpen):
            with patch("os.listdir", new=fake_listdir):
                # noinspection SpellCheckingInspection,SpellCheckingInspection
                class NewFakeFile:
                    """
                    FakeFile class that should fail material lookup
                    """

                    @staticmethod
                    def readline(*_):
                        """
                        readline function that should fail material
                        lookup for thing1
                        """
                        return "ASDF,thing0\n"

                setattr(self.FakeOpen, "FakeFile", NewFakeFile)
                error_string = (
                    "\nMaterial thing1 not found in " + os.path.join(file_directory, "asdfstressweldedasdg") + "\n"
                )
                try:
                    tube._check_materials_list()
                except ValueError as err:
                    tube.TUBE_MATERIALS = old_materials
                    assert str(err) == error_string


def test_import_thread_specs():
    # good input
    test_size = "0-80"
    test_type = ["external", "internal"]
    test_property = "pitch diameter max"
    test_classes = [["2A", "3A"], ["2B", "3B"]]
    good_result = [[0.0514, 0.0519], [0.0542, 0.0536]]
    for thread, classes, expected in zip(test_type, test_classes, good_result):
        # check for correct output
        current_frame = tube.THREAD_SPECS[thread]
        for thread_class, result in zip(classes, expected):
            test_result = current_frame[test_property][test_size][thread_class]
            assert abs(test_result - result) < 1e-7
