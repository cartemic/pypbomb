# -*- coding: utf-8 -*-
"""
Tools for designing and determining operational parameters of a closed-end
detonation tube with optical access.
"""

import os
import warnings
from math import sqrt

import cantera as ct
import numpy as np
import pandas as pd
import pint

from . import thermochem
from . import units
from .thermochem import _U


_Q = _U.Quantity


_DIR_LOOKUP_DATA = os.path.join(
        os.path.dirname(os.path.relpath(__file__)),
        "lookup_data"
)
MATERIAL_LIMITS = dict(
    welded=pd.read_csv(
        os.path.join(
            _DIR_LOOKUP_DATA,
            "ASME_B31_1_stress_limits_welded.csv"
        ),
        # index_col=0
    ),
    seamless=pd.read_csv(
        os.path.join(
            _DIR_LOOKUP_DATA,
            "ASME_B31_1_stress_limits_seamless.csv"
        ),
        # index_col=0
    )
)
PIPE_SCHEDULES = pd.read_csv(
    os.path.relpath(
        os.path.join(
            _DIR_LOOKUP_DATA,
            "pipe_schedules.csv"
        )
    ),
    index_col=0
)


def _collect_tube_materials():
    """
    Reads in a csv file containing tube materials, their corresponding
    ASME B16.5 material groups, and selected material properties.

    Returns
    -------
    df_materials : pd.DataFrame
        Dataframe of materials and their corresponding material groups and
        properties
    """
    file_directory = os.path.join(
        os.path.dirname(
            os.path.relpath(__file__)
        ),
        "lookup_data"
    )
    file_name = "materials_list.csv"
    file_location = os.path.relpath(
        os.path.join(
            file_directory,
            file_name
        )
    )

    # read in csv and extract information
    if os.path.exists(file_location):
        # noinspection PyUnresolvedReferences
        try:
            df_materials = pd.read_csv(file_location)
            # type: pd.DataFrame
        except pd.errors.EmptyDataError:
            raise ValueError("\n" + file_name + " is empty")

    else:
        # raise an exception if the file doesn't exist
        raise ValueError("\n" + file_name + " does not exist")

    # apply units
    df_materials["ElasticModulus"] = [
        _Q(item, "GPa") for item in
        df_materials["ElasticModulus"].values
    ]
    df_materials["Density"] = [
        _Q(item, "g/cm^3") for item in
        df_materials["Density"].values
    ]

    return df_materials


TUBE_MATERIALS = _collect_tube_materials()


def _check_materials_list():
    """
    Makes sure that the materials in materials_list.csv have stress limits
    and flange ratings. This function relies on _get_material_groups(), and
    either raises an error or returns True.

    Returns
    -------
    True
    """
    # collect files
    file_directory = os.path.join(
        os.path.dirname(
            os.path.relpath(__file__)
        ),
        "lookup_data"
    )
    my_files = os.listdir(file_directory)
    flange_ratings = [file for file in my_files if "flange" in file.lower()]
    stress_limits = [file for file in my_files if "stress" in file.lower()]

    # make sure things were actually loaded
    if not bool(flange_ratings + stress_limits):
        raise FileNotFoundError(
            "No files containing 'flange' or 'stress' found"
        )

    # initialize an error string and error indicator. Error string will be
    # used to aggregate errors in the list of available materials so that
    # all issues may be rectified simultaneously.
    error_string = "\n"
    has_errors = False

    # make sure all pipe material limits are either welded or seamless
    # other types are permitted, but will raise a warning
    for file in stress_limits:
        if ("welded" not in file.lower()) and ("seamless" not in file.lower()):
            # warn that something is weird
            warnings.warn(
                file + " does not indicate whether it is welded or seamless"
            )

        # check the first row of the file in question to extract the names
        # of the materials that it contains stress limits for
        file_location = os.path.join(
            file_directory,
            file
        )
        with open(file_location, "r") as current_file:
            # read the first line, strip off carriage return, and split by
            # comma separators. Ignore first value, as this is temperature.
            materials = current_file.readline().strip().split(",")[1:]

            # check to make sure that each material in the list of available
            # materials has a stress limit curve for the current limit type
            for item in TUBE_MATERIALS["Grade"].values.astype(str):
                if item not in materials:
                    # a material is missing from the limits spreadsheet.
                    # indicate that an error has occurred, and add it to the
                    # error string.
                    error_string += "Material " + item + " not found in " \
                                    + file_location + "\n"
                    has_errors = True

    # find out which material groups need to be inspected
    groups = set()
    for group in TUBE_MATERIALS["Group"].values.astype(str):
        groups.add(group.replace(".", "_"))

    # check folder to make sure the correct files exist
    for group in groups:
        if not any(rating.find(group) > 0 for rating in flange_ratings):
            # current group was not found in any of the files
            error_string += "material group " + group + " not found" + "\n"
            has_errors = True

    # report all errors
    if has_errors:
        raise ValueError(error_string)

    return True


_check_materials_list()


def _check_material(material):
    if material not in TUBE_MATERIALS.Grade.values:
        msg = "%s not found in materials." \
              " Try Flange.available_materials" % material
        raise ValueError(msg)
    return True


def _get_flange_limits_from_csv():
    """
    Reads in flange pressure limits as a function of temperature for
    different pressure classes per ASME B16.5. Temperature is in Centigrade
    and pressure is in bar.
    """
    groups = ["2.1", "2.2", "2.3"]
    group_flange_limits = {group: None for group in groups}

    for group in groups:
        # ensure group is valid
        file_group = str(group).replace(".", "_")
        file_directory = os.path.join(
            os.path.dirname(os.path.relpath(__file__)),
            "lookup_data")
        file_name = "ASME_B16_5_flange_ratings_group_" + file_group + \
                    ".csv"
        file_location = os.path.relpath(
            os.path.join(file_directory, file_name)
        )
        if not os.path.exists(file_location):
            raise FileNotFoundError(
                "\n" + file_location + "not found"
            )

        # import the correct .csv file as a pandas dataframe
        flange_limits = pd.read_csv(file_location)

        # ensure all temperatures and pressures are floats
        new_data = np.array([
            pd.to_numeric(flange_limits[column].values, errors="coerce")
            for column in flange_limits.columns
        ]).transpose()
        flange_limits = pd.DataFrame(
            columns=flange_limits.columns,
            data=new_data
        ).fillna(0)

        # make sure pressures are positive
        if not all(
                flange_limits.loc[
                :,
                flange_limits.columns != "Temperature"
                ].fillna(0).values.flatten() >= 0
        ):
            raise ValueError("\nPressure less than zero.")

        # add units to temperature column
        flange_limits["Temperature"] = [
            _Q(temp, "degC") for temp in
            flange_limits["Temperature"]
        ]

        # add units to pressure columns
        for key in flange_limits.keys():
            if key != "Temperature":
                pressures = []
                for pressure in flange_limits[key]:
                    if pressure < 0:
                        pressures.append(np.NaN)
                    else:
                        pressures.append(_Q(
                            float(pressure), "bar")
                        )
                flange_limits[key] = pressures

        group_flange_limits[group] = flange_limits
    return group_flange_limits


FLANGE_LIMITS = _get_flange_limits_from_csv()


def _import_thread_specs():
    """
    Imports thread specifications from .csv files

    Returns
    -------
    thread_specs : dict
        [internal thread specs, external thread specs]. Both sets of thread
        specifications are multi-indexed with (thread size, thread class).
    """
    file_directory = os.path.join(
        os.path.dirname(
            os.path.relpath(__file__)
        ),
        "lookup_data"
    )
    file_names = [
        "ANSI_inch_internal_thread.csv",
        "ANSI_inch_external_thread.csv"
    ]
    file_locations = [
        os.path.relpath(
            os.path.join(
                file_directory,
                name
            )
        )
        for name in file_names
    ]

    thread_specs = {
        key: pd.read_csv(location, index_col=(0, 1)) for location, key in
        zip(file_locations, ["internal", "external"])
    }

    return thread_specs


THREAD_SPECS = _import_thread_specs()


class Bolt:
    """
    Methods relating to bolt calculations and property lookup.
    """
    @classmethod
    def calculate_safety_factors(
            cls,
            max_pressure,
            window_area,
            num_bolts,
            thread_size,
            thread_class,
            bolt_max_tensile,
            plate_max_tensile,
            engagement_length,
            unit_registry=_U
    ):
        """
        Calculates bolt and plate safety factors for viewing window bolts

        Parameters
        ----------
        max_pressure : pint.Quantity or Tuple[float, str]
            Tube maximum pressure
        window_area : pint.Quantity or Tuple[float, str]
            Window area exposed to high pressure environment
        num_bolts : int or np.ndarray
            Number of bolts used to secure each viewing window
        thread_size : str
            Size of threads to be evaluated, e.g. ``1/4-20`` or ``1 1/2-6``
        thread_class : str
            Class of threads to be evaluated, ``"2"`` or ``"3"``. (``"A"`` or
            ``"B"`` are automatically appended for internal/external threads)
        bolt_max_tensile : pint.Quantity or Tuple[float, str]
            Pint quantity of bolt (ext. thread) tensile failure stress
        plate_max_tensile : pint.Quantity or Tuple[float, str]
            Pint quantity of plate (int. thread) tensile failure stress
        engagement_length : pint.Quantity or Tuple[float, str]
            Pint quantity of total thread engagement length
        unit_registry : pint.UnitRegistry
            Keeps output consistent with parent registry, avoiding conflicts

        Returns
        -------
        dict
            Dictionary giving factors of safety for window bolts and the plate

            that they are screwed into. Keys:

            * ``"bolt"``
            * ``"plate"``
        """
        quant = unit_registry.Quantity

        max_pressure = units.parse_quant_input(
            max_pressure,
            unit_registry
        )
        window_area = units.parse_quant_input(
            window_area,
            unit_registry
        )
        bolt_max_tensile = units.parse_quant_input(
            bolt_max_tensile,
            unit_registry
        )
        plate_max_tensile = units.parse_quant_input(
            plate_max_tensile,
            unit_registry
        )
        engagement_length = units.parse_quant_input(
            engagement_length,
            unit_registry
        )

        units.check_pint_quantity(
            max_pressure,
            "pressure",
            ensure_positive=True
        )
        units.check_pint_quantity(
            window_area,
            "area",
            ensure_positive=True
        )
        units.check_pint_quantity(
            bolt_max_tensile,
            "pressure",
            ensure_positive=True
        )
        units.check_pint_quantity(
            plate_max_tensile,
            "pressure",
            ensure_positive=True
        )
        units.check_pint_quantity(
            engagement_length,
            "length",
            ensure_positive=True
        )

        # get total force per bolt
        window_force = (
                (max_pressure - quant(1, "atm")) * window_area / num_bolts
        )

        # get stress areas
        thread = cls.calculate_stress_areas(
            thread_size,
            thread_class,
            bolt_max_tensile,
            plate_max_tensile,
            engagement_length,
            unit_registry
        )
        screw_area = thread["screw area"]
        screw_area = quant(
            screw_area.magnitude,
            screw_area.units.format_babel()
        )
        plate_area = thread["plate area"]
        plate_area = quant(
            plate_area.magnitude,
            plate_area.units.format_babel()
        )

        # calculate safety factors
        safety_factors = dict()
        safety_factors["bolt"] = (
                bolt_max_tensile / (window_force / screw_area)
        ).to("").magnitude
        safety_factors["plate"] = (
                plate_max_tensile / (window_force / plate_area)
        ).to("").magnitude
        return safety_factors

    @staticmethod
    def calculate_stress_areas(
            thread_size,
            thread_class,
            bolt_max_tensile,
            plate_max_tensile,
            engagement_length,
            unit_registry
    ):
        """
        Calculates internal and external thread stress areas using formulas in
        Machinery's Handbook, 26th edition.

        Parameters
        ----------
        thread_size : str
            Size of threads to be evaluated, e.g. ``"1/4-20"`` or ``"1 1/2-6"``
        thread_class : str
            Class of threads to be evaluated, ``"2"`` or ``"3"``. (``"A"`` or
            ``"B"`` are automatically appended for internal/external threads)
        bolt_max_tensile : pint.Quantity or Tuple[float, str]
            Pint quantity of bolt (external thread) tensile failure stress
        plate_max_tensile : pint.Quantity or Tuple[float, str]
            Pint quantity of plate (internal thread) tensile failure stress
        engagement_length : pint.Quantity or Tuple[float, str]
            Pint quantity of total thread engagement length
        unit_registry : pint.UnitRegistry
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        dict
            Dictionary with the following keys:

            * ``"plate area"``: Stress area of internal threads within the plate
            * ``"screw area"``: Stress area of external threads on the screw
            * ``"minimum engagement"``: Minimum engagement length causing screw
              to fail in tension rather than shear, thus preventing the plate
              from stripping.
        """
        quant = unit_registry.Quantity

        units.check_pint_quantity(
            bolt_max_tensile,
            "pressure",
            ensure_positive=True
        )
        units.check_pint_quantity(
            plate_max_tensile,
            "pressure",
            ensure_positive=True
        )
        units.check_pint_quantity(
            engagement_length,
            "length",
            ensure_positive=True
        )

        # convert to local unit registry
        bolt_max_tensile = quant(
            bolt_max_tensile.magnitude,
            bolt_max_tensile.units.format_babel()
        )
        plate_max_tensile = quant(
            plate_max_tensile.magnitude,
            plate_max_tensile.units.format_babel()
        )
        engagement_length = quant(
            engagement_length.magnitude,
            engagement_length.units.format_babel()
        )

        thread = dict()

        # look up thread specs for stress area calculations
        k_n_max = quant(
            THREAD_SPECS["internal"]
            ["minor diameter max"]
            [thread_size]
            [thread_class + "B"],
            "in"
        )
        e_s_min = quant(
            THREAD_SPECS["external"]
            ["pitch diameter min"]
            [thread_size]
            [thread_class + "A"],
            "in"
        )
        e_n_max = quant(
            THREAD_SPECS["internal"]
            ["pitch diameter max"]
            [thread_size]
            [thread_class + "B"],
            "in"
        )
        d_s_min = quant(
            THREAD_SPECS["external"]
            ["major diameter min"]
            [thread_size]
            [thread_class + "A"],
            "in"
        )
        tpi = quant(
            float(thread_size.split("-")[-1]),
            "1/in"
        )
        basic_diameter = quant(
            THREAD_SPECS["external"]
            ["basic diameter"]
            [thread_size]
            [thread_class + "A"],
            "in"
        )

        if bolt_max_tensile < quant(100000, "psi"):
            # calculate screw tensile area using eq. 9 (p. 1482) in Fasteners
            # section of Machinery's Handbook 26 (also eq. 2a on p. 1490)
            screw_area_tensile = np.pi / 4 * (
                    basic_diameter - 0.9742785 / tpi
            ) ** 2
        else:
            # calculate screw tensile area using eq. 2b (p. 1490) in Fasteners
            # section of Machinery's Handbook 26
            screw_area_tensile = np.pi * (
                    e_s_min / 2 -
                    0.16238 / tpi
            ) ** 2

        # calculate screw shear area using eq. 5 (p. 1491) in Fasteners section
        # of Machinery's Handbook 26
        screw_area_shear = (
                np.pi * tpi * engagement_length * k_n_max *
                (1. / (2 * tpi) + 0.57735 * (e_s_min - k_n_max))
        )

        # choose correct area
        if screw_area_shear < screw_area_tensile:
            warnings.warn(
                "Screws fail in shear, not tension." +
                " Plate may be damaged." +
                " Consider increasing bolt engagement length",
                Warning
            )
            thread["screw area"] = screw_area_shear
        else:
            thread["screw area"] = screw_area_tensile

        # calculate plate shear area using eq. 6 (p. 1491) in Fasteners section
        # of Machinery's Handbook 26
        thread["plate area"] = (
                np.pi * tpi * engagement_length * d_s_min *
                (1. / (2 * tpi) + 0.57735 * (d_s_min - e_n_max))
        )

        # calculate minimum engagement scale factor using eq. 3 (p. 1490) in
        # Fasteners section of Machinery's Handbook 26
        j_factor = (
                (screw_area_shear * bolt_max_tensile) /
                (thread["plate area"] * plate_max_tensile)
        )

        # calculate minimum thread engagement (corrected for material
        # differences) using eqs. 1 and 4 (pp. 1490-1491) in Fasteners section
        # of Machinery's Handbook 26
        thread["minimum engagement"] = (
            2 * screw_area_tensile / (
                k_n_max * np.pi * (
                    1. / 2 + 0.57735 * tpi * (e_s_min - k_n_max)
                    )
                )
            ) * j_factor

        return thread

    @staticmethod
    def get_thread_property(
            thread_property,
            thread_size,
            thread_class,
            unit_registry
    ):
        """
        Finds a thread property, such as minor diameter, using a dataframe from
        import_thread_specs(). import_thread_specs is not directly called here
        to save time by not reading from disk every time a property is
        requested.

        Parameters
        ----------
        thread_property : str
            Property that is desired, such as ``"minor diameter"``
        thread_size : str
            Thread size for desired property, such as ``"1/4-20"`` or
            ``"1 1/2-6"``
        thread_class : str
            Thread class: ``"2B"`` or ``"3B"`` for internal threads, ``"2A"`` or
            ``"3A"`` for external threads
        unit_registry : pint.UnitRegistry
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        pint.Quantity
            Property requested, as a pint quantity with units of inches
        """
        quant = unit_registry.Quantity

        # determine if internal or external
        if "A" in thread_class and ("2" in thread_class or "3" in thread_class):
            thread_specs = THREAD_SPECS["external"]
        elif "B" in thread_class and ("2" in thread_class
                                      or "3" in thread_class):
            thread_specs = THREAD_SPECS["internal"]
        else:
            raise ValueError("\nbad thread class")

        # ensure property is in the specs dataframe
        if thread_property not in thread_specs.keys():
            raise KeyError("\nThread property \'" +
                           thread_property +
                           "\' not found. Available specs: " +
                           "'" + "', '".join(thread_specs.keys()) + "'")

        # ensure thread size is in the specs dataframe
        if thread_size not in thread_specs.index:
            raise KeyError("\nThread size \'" +
                           thread_size +
                           "\' not found")

        # retrieve the property
        return quant(thread_specs[thread_property][thread_size][thread_class],
                     "in")


class DDT:
    """
    Methods for estimating the deflagration-to-detonation transition (DDT).
    """
    @staticmethod
    def calculate_blockage_diameter(
            pipe_id,
            blockage_ratio,
            unit_registry=_U
    ):
        """
        Calculates the diameter of a Shchelkin spiral corresponding to a given
        blockage ratio within a pipe of given inner diameter.

        Parameters
        ----------
        pipe_id : pint.Quantity or Tuple[float, str]
            Length scale representing the inner diameter of the pipe used for
            the detonation tube
        blockage_ratio : float
            Ratio of blocked area to total cross-sectional area, :math:`0 < BR
            < 1`
        unit_registry : pint.UnitRegistry
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        pint.Quantity
            Shchelkin spiral diameter inside a tube of `pipe_id` inner diameter
            giving a blockage ratio of `blockage_ratio`. Units are the same as
            pipe_id.
        """
        # ensure blockage ratio is a float
        try:
            blockage_ratio = float(blockage_ratio)
        except ValueError:
            raise ValueError("\nNon-numeric blockage ratio.")

        # ensure blockage ratio is on 0<BR<1
        if not 0 < blockage_ratio < 1:
            raise ValueError("\nBlockage ratio outside of 0<BR<1")

        pipe_id = units.parse_quant_input(pipe_id, unit_registry)
        units.check_pint_quantity(
            pipe_id,
            "length",
            ensure_positive=True
        )

        # calculate Shchelkin spiral diameter
        spiral_diameter = pipe_id / 2 * (1 - sqrt(1 - blockage_ratio))
        return spiral_diameter

    @staticmethod
    def calculate_blockage_ratio(
            tube_inner_diameter,
            blockage_diameter,
            unit_registry=_U
    ):
        """
        Calculates the blockage ratio of a Shchelkin spiral within a detonation
        tube.

        Parameters
        ----------
        tube_inner_diameter : pint.Quantity or Tuple[float, str]
            Inner diameter of the detonation tube
        blockage_diameter : pint.Quantity or Tuple[float, str]
            Outer diameter of the blockage used to create the Shchelkin spiral
            (i.e. a Shchelkin spiral made from 1/2" round stock would be
            ``blockage_diameter=quant(0.5, "inch")``
        unit_registry : pint.UnitRegistry
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        float
            Ratio of blocked to open area (between 0 and 1)
        """

        tube_inner_diameter = units.parse_quant_input(
            tube_inner_diameter,
            unit_registry
        ).to_base_units()
        blockage_diameter = units.parse_quant_input(
            blockage_diameter,
            unit_registry
        ).to_base_units()
        units.check_pint_quantity(
            tube_inner_diameter,
            "length",
            ensure_positive=True
        )
        units.check_pint_quantity(
            blockage_diameter,
            "length",
            ensure_positive=True
        )

        # ensure blockage diameter < tube diameter
        if tube_inner_diameter.magnitude == 0:
            raise ValueError("\ntube ID cannot be 0")
        elif blockage_diameter >= tube_inner_diameter:
            raise ValueError("\nblockage diameter >= tube diameter")

        # calculate blockage ratio
        blockage_ratio = (1 - (1 - 2 * blockage_diameter.magnitude /
                               tube_inner_diameter.magnitude) ** 2)

        return blockage_ratio

    @staticmethod
    def calculate_run_up(
            blockage_ratio,
            tube_diameter,
            initial_temperature,
            initial_pressure,
            species_dict,
            mechanism,
            unit_registry=_U,
            phase_specification="",
    ):
        """
        Calculates the runup distance needed for a detonation to develop from a
        deflagration for a given blockage ratio, tube diameter, and mixture.
        This is accomplished using equations collected by Ciccarelli and
        Dorofeev [1] for blockage ratios <= 0.75. If the desired blockage ratio
        is less than 0.3, the mixture viscosity is needed, and the
        `phase_specification` option may be necessary depending on the
        mechanism.

        [1] G. Ciccarelli and S. Dorofeev, *Flame acceleration and transition to
        detonation in ducts*, Progress in Energy and Combustion Science,
        vol. 34, no. 4, pp. 499–550, Aug. 2008.
        https://doi.org/10.1016/j.pecs.2007.11.002

        Parameters
        ----------
        blockage_ratio : float
            Ratio of blocked area to total cross-sectional area, :math:`0 < BR
            < 1`
        tube_diameter : pint.Quantity or Tuple[float, str]
            Internal diameter of the detonation tube
        initial_temperature : pint.Quantity or Tuple[float, str]
            Mixture initial temperature
        initial_pressure : pint.Quantity or Tuple[float, str]
            Mixture initial pressure
        species_dict : dict
            Dictionary containing the species in the mixture as keys, with total
            moles or mole fractions as values
        mechanism : str
            Mechanism file name for Cantera. See ``tools.find_mechanisms()`` for
            a list of installed mechanisms.
        unit_registry : pint.UnitRegistry
            Unit registry for managing units to prevent conflicts with parent
            unit registry
        phase_specification : str, optional
            Phase specification within the mechanism file used to evaluate
            thermophysical properties. If ``gri30.cti`` is used with no phase
            specification viscosity calculations will fail, resulting in
            an error for all blockage ratios less than 0.3.

        Returns
        -------
        pint.Quantity
            Predicted DDT distance
        """

        if blockage_ratio <= 0 or blockage_ratio > 0.75:
            raise ValueError("\nBlockage ratio outside of correlation range")

        tube_diameter = units.parse_quant_input(
            tube_diameter,
            unit_registry
        )
        initial_temperature = units.parse_quant_input(
            initial_temperature,
            unit_registry
        )
        initial_pressure = units.parse_quant_input(
            initial_pressure,
            unit_registry
        )
        units.check_pint_quantity(
            tube_diameter,
            "length",
            ensure_positive=True
        )
        units.check_pint_quantity(
            initial_temperature,
            "temperature",
            ensure_positive=True
        )
        units.check_pint_quantity(
            initial_pressure,
            "pressure",
            ensure_positive=True
        )

        # calculate laminar flame speed
        laminar_fs = thermochem.calculate_laminar_flame_speed(
            initial_temperature,
            initial_pressure,
            species_dict,
            mechanism,
            phase_specification,
            unit_registry
        )

        # calculate density ratio across the deflagration assuming adiabatic
        # flame
        density = np.zeros(2)
        working_gas = ct.Solution(mechanism, phase_specification)
        working_gas.TPX = [
            initial_temperature.to("K").magnitude,
            initial_pressure.to("Pa").magnitude,
            species_dict
        ]
        density[0] = working_gas.density
        working_gas.equilibrate("HP")
        density[1] = 1 / working_gas.density
        density_ratio = np.prod(density)

        # find sound speed in products at adiabatic flame temperature
        sound_speed = thermochem.get_eq_sound_speed(
            (working_gas.T, "K"),
            (working_gas.P, "Pa"),
            species_dict,
            mechanism,
            phase_specification,
            unit_registry
        )

        def eq4_1():
            """
            Calculate runup distance for blockage ratios <= 0.1 using equation
            4.1 from G. Ciccarelli and S. Dorofeev, “Flame acceleration and
            transition to detonation in ducts,” Progress in Energy and
            Combustion Science, vol. 34, no. 4, pp. 499–550, Aug. 2008.
            """
            # define constants
            kappa = 0.4
            kk = 5.5
            cc = 0.2
            mm = -0.18
            eta = 2.1

            # calculate laminar flame thickness, delta
            working_gas.TPX = [
                initial_temperature.to("K").magnitude,
                initial_pressure.to("Pa").magnitude,
                species_dict
            ]
            rho = units.parse_quant_input(
                (working_gas.density_mass, "kg/m^3"),
                unit_registry
            )
            mu = units.parse_quant_input(
                (working_gas.viscosity, "Pa*s"),
                unit_registry
            )
            nu = mu / rho
            delta = (nu / laminar_fs).to_base_units()

            # calculate gamma
            gamma = (
                            sound_speed /
                            (eta * (density_ratio - 1) ** 2 * laminar_fs) *
                            (delta / tube_diameter) ** (1. / 3)
                    ) ** (1 / (2 * mm + 7. / 3))

            # calculate runup distance
            d_over_h = (2. / (1 - np.sqrt(1 - blockage_ratio)))
            runup = (
                    gamma / cc *
                    (1 / kappa * np.log(gamma * d_over_h) + kk) *
                    tube_diameter
            )
            return runup.to(tube_diameter.units.format_babel())

        def eq4_4():
            """
            Calculate runup for blockage ratios between 0.3 and 0.75 using
            equation 4.4 in G. Ciccarelli and S. Dorofeev, “Flame acceleration
            and transition to detonation in ducts,” Progress in Energy and
            Combustion Science, vol. 34, no. 4, pp. 499–550, Aug. 2008.
            """
            # define constants
            aa = 2.
            bb = 1.5

            # calculate left and right hand sides of eq 4.4
            lhs = (
                    2 * 10 * laminar_fs * (density_ratio - 1) /
                    (sound_speed * tube_diameter)
            )
            rhs = (
                    aa * (1 - blockage_ratio) /
                    (1 + bb * blockage_ratio)
            )

            runup = rhs / lhs

            return runup.to(tube_diameter.units.format_babel())

        # use appropriate equation to calculate runup distance
        if 0.3 <= blockage_ratio <= 0.75:
            runup_distance = eq4_4()
        elif 0.1 >= blockage_ratio:
            runup_distance = eq4_1()
        else:
            interp_distances = np.array([
                eq4_1().magnitude,
                eq4_4().magnitude
            ])
            runup_distance = np.interp(
                blockage_ratio,
                np.array([0.1, 0.3]),
                interp_distances
            )
            runup_distance = units.parse_quant_input(
                (runup_distance, tube_diameter.units.format_babel()),
                unit_registry
            )

        return runup_distance


class Window:
    """
    Methods for the design of viewing windows for optical access.
    """
    @classmethod
    def safety_factor(
            cls,
            length,
            width,
            thickness,
            pressure,
            rupture_modulus,
            unit_registry=_U
    ):
        """
        Calculates the safety factor of a clamped rectangular
        window given window dimensions, design pressure, and material rupture
        modulus

        Parameters
        ----------
        length : pint.Quantity or Tuple[float, str]
            Window unsupported (viewing) length
        width : pint.Quantity or Tuple[float, str]
            Window unsupported (viewing) width
        thickness : pint.Quantity or Tuple[float, str]
            Window thickness
        pressure : pint.Quantity or Tuple[float, str]
            Design pressure differential across window at which factor of
            safety is to be calculated
        rupture_modulus : pint.Quantity or Tuple[float, str]
            Rupture modulus of desired window material
        unit_registry : pint.UnitRegistry
            Keeps output consistent with parent registry, avoiding conflicts

        Returns
        -------
        float
            Window factor of safety
        """
        length = units.parse_quant_input(
            length,
            unit_registry
        ).to_base_units()
        width = units.parse_quant_input(
            width,
            unit_registry
        ).to_base_units()
        thickness = units.parse_quant_input(
            thickness,
            unit_registry
        ).to_base_units()
        pressure = units.parse_quant_input(
            pressure,
            unit_registry
        ).to_base_units()
        rupture_modulus = units.parse_quant_input(
            rupture_modulus,
            unit_registry
        ).to_base_units()

        units.check_pint_quantity(
            length,
            "length",
            ensure_positive=True
        )
        units.check_pint_quantity(
            width,
            "length",
            ensure_positive=True
        )
        units.check_pint_quantity(
            thickness,
            "length",
            ensure_positive=True
        )

        units.check_pint_quantity(
            pressure,
            "pressure",
            ensure_positive=True
        )
        units.check_pint_quantity(
            rupture_modulus,
            "pressure",
            ensure_positive=True
        )

        # https://www.crystran.co.uk/userfiles/files/design-of-pressure-windows.pdf
        k = 0.75
        r = length / width
        safety_factor = (
                (thickness / length) ** 2 *
                (2 * rupture_modulus * (1 + r**2)) /
                (k * pressure)
        ).to("").magnitude
        return safety_factor

    @classmethod
    def minimum_thickness(
            cls,
            length,
            width,
            safety_factor,
            pressure,
            rupture_modulus,
            unit_registry=_U
    ):
        """
        Calculates the thickness of a clamped rectangular window which gives
        the desired safety factor

        Parameters
        ----------
        length : pint.Quantity or Tuple[float, str]
            Window unsupported (viewing) length
        width : pint.Quantity or Tuple[float, str]
            Window unsupported (viewing) width
        safety_factor : float
            Safety factor
        pressure : pint.Quantity or Tuple[float, str]
            Design pressure differential across window at which factor of
            safety is to be calculated
        rupture_modulus : pint.Quantity or Tuple[float, str]
            Rupture modulus of desired window material
        unit_registry : pint.UnitRegistry
            Keeps output consistent with parent registry, avoiding conflicts

        Returns
        -------
        pint.Quantity
            Window thickness
        """
        length = units.parse_quant_input(
            length,
            unit_registry
        ).to_base_units()
        width = units.parse_quant_input(
            width,
            unit_registry
        ).to_base_units()
        pressure = units.parse_quant_input(
            pressure,
            unit_registry
        ).to_base_units()
        rupture_modulus = units.parse_quant_input(
            rupture_modulus,
            unit_registry
        ).to_base_units()

        units.check_pint_quantity(
            length,
            "length",
            ensure_positive=True
        )
        units.check_pint_quantity(
            width,
            "length",
            ensure_positive=True
        )
        units.check_pint_quantity(
            pressure,
            "pressure",
            ensure_positive=True
        )
        units.check_pint_quantity(
            rupture_modulus,
            "pressure",
            ensure_positive=True
        )

        # Ensure safety factor is numeric and > 1
        try:
            if safety_factor < 1:
                raise ValueError("\nWindow safety factor < 1")
        except TypeError:
            raise TypeError("\nNon-numeric window safety factor")

        # https://www.crystran.co.uk/userfiles/files/design-of-pressure-windows.pdf
        k = 0.75  # clamped
        r = length / width
        thickness = (
                length *
                np.sqrt(safety_factor * k / 2) *
                np.sqrt(pressure / (rupture_modulus * (1 + r**2)))
        ).to(width.units.format_babel())

        return thickness


class Tube:
    """
    Methods for designing and determining operational limits of a detonation
    tube.
    """
    available_pipe_sizes = list(PIPE_SCHEDULES.index.values)
    available_materials = list(TUBE_MATERIALS.Grade.values)

    @classmethod
    def calculate_max_stress(
            cls,
            initial_temperature,
            material,
            welded,
            unit_registry=thermochem._U
    ):
        """
        Finds the maximum allowable stress of a tube material at the tube's
        initial temperature

        Parameters
        ----------
        initial_temperature : pint.Quantity or Tuple[float, str]
        material : str
            Pipe material (see ``Tube.available_materials``)
        welded : bool
            ``True`` for welded pipe; ``False`` for seamless
        unit_registry : pint.UnitRegistry
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        pint.Quantity
            Maximum allowable tube stress
        """
        initial_temperature = units.parse_quant_input(
            initial_temperature,
            unit_registry
        )

        # look up stress-temperature limits and units
        stress_limits = cls._get_pipe_stress_limits(
            material=material,
            welded=welded,
            unit_registry=unit_registry
        )
        stress_units = stress_limits.values[0].units
        stresses = list(
            map(
                lambda x: x.magnitude,
                stress_limits.values
            )
        )
        temp_units = stress_limits.index[0].units
        temperatures = list(
            map(
                lambda x: x.magnitude,
                stress_limits.index.values
            )
        )

        # ensure material stress limits have monotonically increasing
        # temperatures, otherwise the np.interp "results are nonsense" per
        # scipy docs
        if not np.all(np.diff(temperatures) > 0):
            raise ValueError("\nStress limits require temperatures to be " +
                             "monotonically increasing")

        # interpolate max stress
        max_stress = unit_registry.Quantity(
            np.interp(
                initial_temperature.to(temp_units).magnitude,
                temperatures,
                stresses
            ),
            stress_units
        )

        return max_stress

    @staticmethod
    def _get_pipe_stress_limits(
            material,
            welded,
            unit_registry=_U
    ):
        """
        Looks up ASME B31.1 stress limits as a function of pipe temperature

        Parameters
        ----------
        material : str
            Pipe material (see ``Tube.available_materials``)
        welded : bool
            ``True`` for welded pipe; ``False`` for seamless
        unit_registry : pint.UnitRegistry, optional
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        pd.Series
            Pandas series with data consisting of tress limits (as pint
            quantities), and indices of the corresponding temperatures (also
            as pint quantities)
        """
        quant = unit_registry.Quantity
        if welded:
            material_limits = MATERIAL_LIMITS["welded"][["Temp", material]]
        else:
            material_limits = MATERIAL_LIMITS["seamless"][["Temp", material]]
        material_limits.columns = ["temperature", "stress"]

        material_limits = pd.Series(
            data=[quant(t, "ksi") for t in material_limits["stress"].values],
            index=[quant(t, "degF")
                   for t in material_limits["temperature"].values],
            name="stress"
        )

        return material_limits

    @staticmethod
    def dynamic_load_factor(
            tube_id,
            tube_od,
            cj_velocity,
            elastic_modulus,
            density,
            poisson_ratio,
            plus_or_minus=0.1
    ):
        """
        Calculates the dynamic load factor (DLF) by which a detonation tube's
        static analysis should be scaled in order to account for the tube's
        response to transient pressures. DLF calculation is based on the work of
        Shepherd [1]. Since the limits of "approximately equal to" are not
        defined in the paper, a default value of plus or minus ten percent
        is assumed, thus `plus_or_minus=0.1`.

        [1] Shepherd, J. E. (2009). *Structural Response of Piping to
        Internal Gas Detonation*. Journal of Pressure Vessel Technology,
        131(3), 031204. https://doi.org/10.1115/1.3089497

        Parameters
        ----------
        tube_id : pint.Quantity or Tuple[float, str]
            Inner diameter of tube
        tube_od : pint.Quantity or Tuple[float, str]
            Outer diameter of tube
        cj_velocity : pint.Quantity or Tuple[float, str]
            CJ velocity
        elastic_modulus : pint.Quantity or Tuple[float, str]
            Elastic modulus of tube material
        density : pint.Quantity or Tuple[float, str]
            Density of tube material
        poisson_ratio : float
            Poisson ratio of tube material
        plus_or_minus : float, optional
            Defines the band about the critical velocity which is considered
            "approximately equal to" -- the default value of 0.1 means plus
            or minus ten percent

        Returns
        -------
        float
            Factor by which the tube's static maximum pressure should be
            de-rated to account for transient response to detonation waves
        """
        units.check_pint_quantity(tube_id, "length", True)
        units.check_pint_quantity(tube_od, "length", True)
        units.check_pint_quantity(cj_velocity, "velocity")
        units.check_pint_quantity(elastic_modulus, "pressure", True)
        units.check_pint_quantity(density, "density", True)

        if not (0 < plus_or_minus < 1):
            raise ValueError(
                "\nplus_or_minus factor not between 0 and 1"
            )

        # set geometry
        # first /2 for averaging
        # second /2 to to convert diameter to radius
        radius = (tube_od + tube_id) / 4.
        tube_thickness = (tube_od - tube_id) / 2.

        # calculate critical velocity
        crit_velocity = (
                ((elastic_modulus ** 2 * tube_thickness
                  ** 2) /
                 (3. * density ** 2 * radius ** 2 * (1. - poisson_ratio ** 2))
                 ) ** (1. / 4)
        )

        # set limits for "approximately Vcrit"
        bounds = [
            (crit_velocity * (1. + plus_or_minus)).to("m/s").magnitude,
            (crit_velocity * (1. - plus_or_minus)).to("m/s").magnitude
        ]
        cj_velocity = cj_velocity.to("m/s").magnitude

        if cj_velocity < bounds[1]:
            dynamic_load_factor = 1
        elif cj_velocity > bounds[0]:
            dynamic_load_factor = 2
        else:
            dynamic_load_factor = 4

        return dynamic_load_factor

    @staticmethod
    def calculate_max_pressure(
            tube_id,
            tube_od,
            max_stress,
            safety_factor=1
    ):
        """
        Calculates the maximum allowable pressure from the tube dimensions
        and stress limits using the basic longitudinal joint formula [1].

        [1] E. F. Megyesy, *Pressure vessel handbook*, Oklahoma City, OK:
        PV Publishing, Inc., 2001, p. 14.

        Parameters
        ----------
        tube_id : pint.Quantity or Tuple[float, str]
            Inner diameter of tube
        tube_od : pint.Quantity or Tuple[float, str]
            Outer diameter of tube
        max_stress : pint.Quantity or Tuple[float, str]
            Maximum allowable stress in tube material (e.g. from ASME B31.1)
        safety_factor : float, optional
            Desired safety factor

        Returns
        -------
        pint.Quantity
            Pressure resulting in maximum allowable stress
        """
        mean_diameter = (tube_od + tube_id) / 2.
        wall_thickness = (tube_od - tube_id) / 2.
        max_pressure = max_stress * (2 * wall_thickness) / \
            (mean_diameter * safety_factor)

        return max_pressure

    @classmethod
    def calculate_max_initial_pressure(
            cls,
            tube_id,
            tube_od,
            initial_temperature,
            species_dict,
            mechanism,
            max_pressure,
            elastic_modulus,
            density,
            poisson_ratio,
            plus_or_minus=0.1,
            use_dynamic_load_factor=True,
            unit_registry=pint.UnitRegistry(),
            error_tol=1e-4,
            max_iterations=500,
            use_multiprocessing=False,
            max_pressure_state="reflected",
            return_dlf=False
    ):
        """
        Parameters
        ----------
        tube_id : pint.Quantity or Tuple[float, str]
            Inner diameter of tube
        tube_od : pint.Quantity or Tuple[float, str]
            Outer diameter of tube
            Internal diameter of the detonation tube
        initial_temperature : pint.Quantity or Tuple[float, str]
            Mixture initial temperature
        species_dict : dict
            Dictionary containing the species in the mixture as keys, with total
            moles or mole fractions as values
        mechanism : str
            Mechanism file name for Cantera
        max_pressure : pint.Quantity or Tuple[float, str]
            Maximum allowable pressure within the tube
        elastic_modulus : pint.Quantity or Tuple[float, str]
            Elastic modulus of tube material
        density : pint.Quantity or Tuple[float, str]
            Density of tube material
        poisson_ratio : float
            Poisson ratio of tube material
        plus_or_minus : float
            For dynamic load factor calculation, ``plus_or_minus`` defines
            the band about the critical velocity which is considered
            "approximately equal to." The default value of 0.1 means plus
            or minus ten percent.
        use_dynamic_load_factor : bool
            ``True`` to use dynamic load factor, ``False`` to ignore it.
        unit_registry : pint unit registry
            Unit registry for managing units to prevent conflicts with parent
            unit registry
        error_tol : float
            Relative error tolerance below which initial pressure calculations
            are considered 'good enough'
        max_iterations : int
            Maximum number of loop iterations before exit, defaults to 500
        use_multiprocessing : bool
            ``True`` to use multiprocessing. Defaults to ``False``.
        max_pressure_state : str
            Which state to use as the tube's maximum pressure. Must be either
            ``"reflected"`` (closed tube) or ``"cj"`` (open tube).
        return_dlf : bool
            If ``True``, dynamic load factor will be returned

        Returns
        -------
        pint.Quantity or Tuple[pint.Quantity, float]
            Initial mixture pressure corresponding to the tube's maximum
            allowable pressure.

            If `return_dlf` is set to ``True``, this will be a list, with item
            0 being the initial pressure and item 1 being DLF
        """
        if max_pressure_state not in ("reflected", "cj"):
            raise ValueError("`max_pressure_state` must be 'reflected' or 'cj'")

        quant = unit_registry.Quantity

        # get a rough estimate of the initial pressure
        initial_pressure = quant(1, "atm")
        counter = 0
        error_tol = abs(error_tol)

        state = thermochem.calculate_reflected_shock_state(
            initial_temperature,
            initial_pressure,
            species_dict,
            mechanism,
            unit_registry,
            use_multiprocessing
        )

        dlf = cls.dynamic_load_factor(
            tube_id,
            tube_od,
            state["cj"]["speed"],
            elastic_modulus,
            density,
            poisson_ratio,
            plus_or_minus
        )

        error = (initial_pressure.magnitude * dlf / max_pressure.magnitude) - 1.
        while abs(error) > error_tol and counter < max_iterations:
            counter += 1

            # update initial pressure guess
            initial_pressure = initial_pressure * max_pressure.magnitude / \
                (dlf * state["reflected"]["state"].P)

            # get reflected shock pressure
            state = thermochem.calculate_reflected_shock_state(
                initial_temperature,
                initial_pressure,
                species_dict,
                mechanism,
                unit_registry,
                use_multiprocessing
            )

            if use_dynamic_load_factor:
                dlf = cls.dynamic_load_factor(
                    tube_id,
                    tube_od,
                    state["cj"]["speed"],
                    elastic_modulus,
                    density,
                    poisson_ratio,
                    plus_or_minus
                )
            else:
                dlf = 1.

            error = (state[max_pressure_state]["state"].P * dlf -
                     max_pressure.magnitude) / max_pressure.magnitude

        if return_dlf:
            return initial_pressure, dlf
        else:
            return initial_pressure

    @staticmethod
    def get_available_pipe_schedules(pipe_size):
        """
        Gets available pipe schedules for a given nominal size

        Parameters
        ----------
        pipe_size : str
            Nominal pipe size (see ``Tube.available_pipe_sizes``)

        Returns
        -------
        list
            List of available pipe schedules
        """
        if pipe_size not in PIPE_SCHEDULES.index:
            msg = "Invalid pipe size: %s. " \
                  "See Tube.available_pipe_sizes." % str(pipe_size)
            raise ValueError(msg)
        else:
            return list(PIPE_SCHEDULES.T[pipe_size].dropna().index.drop("OD"))

    @classmethod
    def get_dimensions(
            cls,
            pipe_size,
            pipe_schedule,
            unit_registry=_U
    ):
        """

        Parameters
        ----------
        pipe_size : str
            Nominal pipe size (see ``Tube.available_pipe_sizes``)
        pipe_schedule : str
            Pipe schedule (see ``Tube.get_available_pipe_schedules`` for a list
            of available schedules for `pipe_size`)
        unit_registry : pint.UnitRegistry, optional
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        dict
            Dictionary of pipe dimensions with the keys:

            * ``"inner_diameter"``
            * ``"outer_diameter"``
            * ``"wall_thickness"``
        """
        # note: this also checks for valid size
        if pipe_schedule not in cls.get_available_pipe_schedules(pipe_size):
            msg = "Schedule {:s} invalid for pipe size {:s}. See " \
                  "Tube.available_pipe_sizes and " \
                  "Tube.get_available_pipe_schedules.".format(
                      str(pipe_schedule),
                      str(pipe_size)
                  )
            raise ValueError(msg)
        outer_diameter = PIPE_SCHEDULES["OD"][pipe_size]
        wall_thickness = PIPE_SCHEDULES[pipe_schedule][pipe_size]
        inner_diameter = outer_diameter - 2 * wall_thickness

        return dict(
            inner_diameter=unit_registry.Quantity(inner_diameter, "in"),
            outer_diameter=unit_registry.Quantity(outer_diameter, "in"),
            wall_thickness=unit_registry.Quantity(wall_thickness, "in"),
        )

    @classmethod
    def get_elastic_modulus(
            cls,
            material,
            unit_registry=_U
    ):
        """
        Gets the elastic modulus for a valid tube material. For valid materials
        see ``Tube.available_materials``.

        Parameters
        ----------
        material : str
            Pipe material (see ``Tube.available_materials``)
        unit_registry : pint.UnitRegistry, optional
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        pint.Quantity
            Elastic modulus
        """
        cls._check_material(material)
        return units.parse_quant_input(
            TUBE_MATERIALS[
                TUBE_MATERIALS["Grade"] == material
            ]["ElasticModulus"].values[0],
            unit_registry
        )

    @classmethod
    def get_density(
            cls,
            material,
            unit_registry=_U
    ):
        """
        Gets the density for a valid tube material. For valid materials see
        ``Tube.available_materials``.

        Parameters
        ----------
        material : str
            Pipe material (see ``Tube.available_materials``)
        unit_registry : pint.UnitRegistry, optional
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        pint.Quantity
            Density
        """
        cls._check_material(material)
        return units.parse_quant_input(
            TUBE_MATERIALS[
                TUBE_MATERIALS["Grade"] == material
            ]["Density"].values[0],
            unit_registry
        )

    @classmethod
    def get_poisson(
            cls,
            material
    ):
        """
        Gets the Poisson ratio for a valid tube material. For valid materials
        see ``Tube.available_materials``.

        Parameters
        ----------
        material : str
            Pipe material (see ``Tube.available_materials``)

        Returns
        -------
        float
            Poisson ratio (unitless)
        """
        cls._check_material(material)
        return TUBE_MATERIALS[
            TUBE_MATERIALS["Grade"] == material
        ]["Poisson"].values[0]

    @classmethod
    def _check_material(
            cls,
            material
    ):
        """
        Checks a material to make sure it is in ``Tube.available_materials``.

        Parameters
        ----------
        material

        Returns
        -------
        """
        material = str(material)
        if material not in cls.available_materials:
            raise ValueError("Material %s not found" % material)


class Flange:
    """
    Methods pertaining to flange class and max pressure calculations
    """
    available_materials = list(TUBE_MATERIALS.Grade[
        pd.notna(TUBE_MATERIALS.Group)
    ].values)

    @staticmethod
    def _check_flange_class(flange_class, group):
        if flange_class not in FLANGE_LIMITS[group].keys():
            msg = "{:s} not found in group {:s}.".format(
                str(flange_class),
                str(group),
            )
            raise ValueError(msg)
        return True

    @staticmethod
    def get_class(
            max_pressure,
            temperature,
            material,
            unit_registry=_U
    ):
        """
        Finds the minimum allowable flange class per ASME B16.5 for a given
        flange temperature and tube pressure.

        Parameters
        ----------
        max_pressure : pint.Quantity or Tuple[float, str]
            Maximum pressure within the system as a quantity or tuple of
            ``(magnitude, "units")``
        temperature : pint.Quantity or Tuple[float, str]
            Flange temperature as a quantity or tuple of
            ``(magnitude, "units")``
        material : str
            Flange material (see ``Flange.available_materials``)
        unit_registry : pint.UnitRegistry, optional
            Pint unit registry, if output within a particular registry is
            desired

        Returns
        -------
        str
            Minimum allowable flange class
        """
        max_pressure = units.parse_quant_input(
            max_pressure,
            unit_registry
        )
        temperature = units.parse_quant_input(
            temperature,
            unit_registry
        )
        units.check_pint_quantity(
            max_pressure,
            "pressure",
            True
        )
        units.check_pint_quantity(
            temperature,
            "temperature",
            True
        )
        max_pressure = units.parse_quant_input(
            max_pressure,
            _U
        )
        temperature = units.parse_quant_input(
            temperature,
            _U
        )

        _check_material(material)
        group = str(
            TUBE_MATERIALS[TUBE_MATERIALS.Grade == material].Group.values[0]
        )
        df_limits = FLANGE_LIMITS[group]

        # locate max pressure
        class_keys = df_limits.keys()[1:]
        max_key = "0"
        for key in class_keys:
            if int(key) > int(max_key):
                max_key = key
        max_ok_pressure = df_limits[max_key].dropna().max()

        # ensure pressure is within bounds
        if max_pressure > max_ok_pressure:
            raise ValueError("Pressure out of range.")

        # ensure temperature is within bounds
        if (
                temperature < df_limits["Temperature"].min()
        ) or (
                temperature > df_limits["Temperature"].max()
        ):
            raise ValueError("Temperature out of range.")

        df_limits = df_limits.applymap(
            lambda x: x.to_base_units().magnitude
        )
        df_interp = pd.DataFrame(columns=df_limits.columns)
        df_interp["Temperature"] = [temperature.to_base_units().magnitude]
        df_interp = pd.concat((df_limits, df_interp)).set_index("Temperature")
        with warnings.catch_warnings():
            # ignore max iteration warning during interpolation
            warnings.simplefilter("ignore")
            df_interp.interpolate(method="slinear", inplace=True)
        ser_result = df_interp.iloc[-1]
        return ser_result.keys()[
            ser_result >= max_pressure.to_base_units().magnitude
        ][0]

    @staticmethod
    def get_max_pressure(
            flange_class,
            temperature,
            material
    ):
        """
        Finds the max allowable pressure for a flange of a given material and
        class at the desired temperature.

        Parameters
        ----------
        flange_class : str
            Class of the flange to evaluate, e.g. ``"900"``
        temperature : pint.Quantity or Tuple[float, str]
            Flange temperature
        material : str
            Flange material (see ``Flange.available_materials``)

        Returns
        -------
        pint.Quantity
            Maximum allowable pressure
        """
        _check_material(material)
        group = str(TUBE_MATERIALS.Group[
            TUBE_MATERIALS.Grade == material
        ].iloc[0])
        Flange._check_flange_class(flange_class, group)

        # ensure correct unit registry is used for comparison
        if hasattr(temperature, "_REGISTRY"):
            ureg_out = temperature._REGISTRY
        else:
            ureg_out = _U

        temperature = units.parse_quant_input(temperature, _U)

        df_limits = FLANGE_LIMITS[group]
        # ensure temperature is in range
        if (
                temperature < df_limits["Temperature"].min()
        ) or (
                temperature > df_limits["Temperature"].max()
        ):
            raise ValueError("Temperature out of range.")

        df_limits = df_limits.applymap(
            lambda x: x.to_base_units().magnitude
        )

        df_interp = pd.DataFrame(columns=df_limits.columns)
        df_interp["Temperature"] = [temperature.to_base_units().magnitude]
        df_interp = pd.concat((df_limits, df_interp)).set_index("Temperature")
        with warnings.catch_warnings():
            # ignore max iteration warning during interpolation
            warnings.simplefilter("ignore")
            df_interp.interpolate(method="slinear", inplace=True)
        ser_result = df_interp.iloc[-1]
        return ureg_out.Quantity(ser_result[flange_class], "Pa")

    @staticmethod
    def get_max_temperature(
            flange_class,
            pressure,
            material
    ):
        """
        Finds the max allowable temperature for a flange of a given material and
        class at the desired pressure.

        Parameters
        ----------
        flange_class : str
            Class of the flange to evaluate, e.g. ``"900"``
        pressure : pint.Quantity or Tuple[float, str]
            System pressure
        material : str
            Flange material (see ``Flange.available_materials``)

        Returns
        -------
        pint.Quantity
            Max allowable flange temperature
        """
        _check_material(material)
        group = str(TUBE_MATERIALS.Group[
            TUBE_MATERIALS.Grade == material
        ].iloc[0])
        Flange._check_flange_class(flange_class, group)

        # ensure correct unit registry is used for comparison
        if hasattr(pressure, "_REGISTRY"):
            ureg_out = pressure._REGISTRY
        else:
            ureg_out = _U

        pressure = units.parse_quant_input(pressure, _U)

        df_limits = FLANGE_LIMITS[group]
        # ensure pressure is in range
        if (
                pressure < df_limits[flange_class].min()
        ) or (
                pressure > df_limits[flange_class].max()
        ):
            raise ValueError("Pressure out of range.")

        df_limits = df_limits.applymap(
            lambda x: x.to_base_units().magnitude
        )[["Temperature", flange_class]]

        df_interp = pd.DataFrame(columns=df_limits.columns)
        df_interp[flange_class] = [pressure.to_base_units().magnitude]
        df_interp = pd.concat((df_limits, df_interp)).set_index(flange_class)
        with warnings.catch_warnings():
            # ignore max iteration warning during interpolation
            warnings.simplefilter("ignore")
            df_interp.interpolate(method="slinear", inplace=True)
        ser_result = df_interp.iloc[-1]
        return ureg_out.Quantity(ser_result["Temperature"], "K")
