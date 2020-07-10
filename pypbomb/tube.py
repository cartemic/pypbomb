# -*- coding: utf-8 -*-
"""
PURPOSE:
    A series of tools to aid in the design of a detonation tube.

CREATED BY:
    Mick Carter
    Oregon State University
    CIRE and Propulsion Lab
    cartemic@oregonstate.edu
"""

import os
import warnings
from math import sqrt

import cantera as ct
import numpy as np
import pandas as pd
import pint
import sympy as sp

from . import thermochem
from . import tools
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


class Bolt:
    @classmethod
    def calculate_stress_areas(
            cls,
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
            Size of threads to be evaluated, e.g. '1/4-20' or '1 1/2-6'
        thread_class : str
            Class of threads to be evaluated, '2' or '3'. 'A' or 'B' are
            automatically appended for internal/external threads
        bolt_max_tensile : pint quantity
            Pint quantity of bolt (ext. thread) tensile failure stress
        plate_max_tensile : pint quantity
            Pint quantity of plate (int. thread) tensile failure stress
        engagement_length : pint quantity
            Pint quantity of total thread engagement length
        unit_registry : pint unit registry
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        thread : dict
            Dictionary with the following key/value pairs:
            'plate area': stress area of internal threads within the plate
            'screw area': stress area of external threads on the screw
            'minimum engagement': minimum engagement length causing screw to
                fail in tension rather than shear, thus preventing the plate
                from stripping.
        """
        quant = unit_registry.Quantity

        tools.check_pint_quantity(
            bolt_max_tensile,
            "pressure",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            plate_max_tensile,
            "pressure",
            ensure_positive=True
        )
        tools.check_pint_quantity(
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
        thread_specs = cls._import_thread_specs()
        k_n_max = quant(
            thread_specs["internal"]
            ["minor diameter max"]
            [thread_size]
            [thread_class + "B"],
            "in"
        )
        e_s_min = quant(
            thread_specs["external"]
            ["pitch diameter min"]
            [thread_size]
            [thread_class + "A"],
            "in"
        )
        e_n_max = quant(
            thread_specs["internal"]
            ["pitch diameter max"]
            [thread_size]
            [thread_class + "B"],
            "in"
        )
        d_s_min = quant(
            thread_specs["external"]
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
            thread_specs["external"]
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

    @classmethod
    def get_thread_property(
            cls,
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
            Property that is desired, such as 'minor diameter'
        thread_size : str
            Thread size for desired property, such as '1/4-20' or '1 1/2-6'
        thread_class : str
            Thread class: '2B' or '3B' for internal threads, '2A' or '3A' for
            external threads
        unit_registry : pint unit registry
            Unit registry for managing units to prevent conflicts with parent
            unit registry

        Returns
        -------
        pint.UnitRegistry().Quantity
            Property requested, as a pint quantity with units of inches
        """
        quant = unit_registry.Quantity
        thread_specs = cls._import_thread_specs()

        # determine if internal or external
        if "A" in thread_class and ("2" in thread_class or "3" in thread_class):
            thread_specs = thread_specs["external"]
        elif "B" in thread_class and ("2" in thread_class
                                      or "3" in thread_class):
            thread_specs = thread_specs["internal"]
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
    @staticmethod
    def calculate_spiral_diameter(
            pipe_id,
            blockage_ratio
    ):
        """
        Calculates the diameter of a Shchelkin spiral corresponding to a given
        blockage ratio within a pipe of given inner diameter.

        Parameters
        ----------
        pipe_id : pint quantity
            Length scale representing the inner diameter of the pipe used for
            the detonation tube
        blockage_ratio : float
            percentage (float between 0 and 1)

        Returns
        -------
        spiral_diameter : pint quantity
            Shchelkin spiral diameter inside a tube of pipe_id inner diameter
            giving a blockage ratio of blockage_ratio %. Units are the same as
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

        tools.check_pint_quantity(
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
            blockage_diameter
    ):
        """
        Calculates the blockage ratio of a Shchelkin spiral within a detonation
        tube.

        Parameters
        ----------
        tube_inner_diameter : pint quantity
            Length scale corresponding to the ID of the detonation tube
        blockage_diameter : pint quantity
            Length scale corresponding to the OD of a Shchelkin spiral

        Returns
        -------
        blockage_ratio : float
            Ratio of blocked to open area (between 0 and 1)
        """

        # check dimensionality and >=0
        tools.check_pint_quantity(
            tube_inner_diameter,
            "length",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            blockage_diameter,
            "length",
            ensure_positive=True
        )

        # make sure units cancel
        blockage_diameter = blockage_diameter.to_base_units()
        tube_inner_diameter = tube_inner_diameter.to_base_units()

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
            unit_registry,
            phase_specification=""
    ):
        """
        Calculates the runup distance needed for a detonation to develop from a
        deflagration for a given blockage ratio, tube diameter, and mixture.
        This is accomplished using equations collected by Ciccarelli and
        Dorofeev [1] for blockage ratios <= 0.75. If the desired blockage ratio
        is less than 0.3, the mixture viscosity is needed, and the
        phase_specification option may be necessary depending on the mechanism.

        [1] G. Ciccarelli and S. Dorofeev, “Flame acceleration and transition to
        detonation in ducts,” Progress in Energy and Combustion Science,
        vol. 34, no. 4, pp. 499–550, Aug. 2008.

        Parameters
        ----------
        blockage_ratio : float
            Ratio of the cross-sectional area of the detonation tube and a
            periodic blockage used to cause DDT
        tube_diameter : pint quantity
            Internal diameter of the detonation tube
        initial_temperature : pint quantity
            Mixture initial temperature
        initial_pressure : pint quantity
            Mixture initial pressure
        species_dict : dict
            Dictionary containing the species in the mixture as keys, with total
            moles or mole fractions as values
        mechanism : str
            Mechanism file name for Cantera
        unit_registry : pint unit registry
            Unit registry for managing units to prevent conflicts with parent
            unit registry
        phase_specification : str
            (Optional) Phase specification within the mechanism file used to
            evaluate thermophysical properties. If Gri30.cti is used with no
            phase specification, viscosity calculations will fail, resulting in
            an error for all blockage ratios less than 0.3.

        Returns
        -------
        runup_distance : pint quantity
            Predicted DDT distance, with the same units as the tube diameter
        """

        if blockage_ratio <= 0 or blockage_ratio > 0.75:
            raise ValueError("\nBlockage ratio outside of correlation range")

        tools.check_pint_quantity(
            tube_diameter,
            "length",
            ensure_positive=True
        )

        tools.check_pint_quantity(
            initial_temperature,
            "temperature",
            ensure_positive=True
        )

        tools.check_pint_quantity(
            initial_pressure,
            "pressure",
            ensure_positive=True
        )

        # handle units
        quant = unit_registry.Quantity
        tube_diameter = quant(
            tube_diameter.magnitude,
            tube_diameter.units.format_babel()
        )

        # calculate laminar flame speed
        laminar_fs = thermochem.calculate_laminar_flame_speed(
            initial_temperature,
            initial_pressure,
            species_dict,
            mechanism
        )
        laminar_fs = quant(
            laminar_fs.magnitude, laminar_fs.units.format_babel()
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
            quant(working_gas.T, "K"),
            quant(working_gas.P, "Pa"),
            species_dict,
            mechanism
        )
        sound_speed = quant(
            sound_speed.magnitude,
            sound_speed.units.format_babel()
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
            rho = quant(working_gas.density_mass, "kg/m^3")
            mu = quant(working_gas.viscosity, "Pa*s")
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
            runup_distance = quant(
                runup_distance,
                tube_diameter.units.format_babel()
            )

        return runup_distance


class Window:
    @classmethod
    def safety_factor(
            cls,
            length,
            width,
            thickness,
            pressure,
            rupture_modulus
    ):
        """
        This function calculates the safety factor of a clamped rectangular
        window given window dimensions, design pressure, and material rupture
        modulus

        Parameters
        ----------
        length : pint quantity with length units
            Window unsupported (viewing) length
        width : pint quantity with length units
            Window unsupported (viewing) width
        thickness : pint quantity with length units
            Window thickness
        pressure : pint quantity with pressure units
            Design pressure differential across window at which factor of
            safety is to be calculated
        rupture_modulus : pint quantity with pressure units
            Rupture modulus of desired window material.

        Returns
        -------
        safety_factor : float
            Window factor of safety
        """

        tools.check_pint_quantity(
            length,
            "length",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            width,
            "length",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            thickness,
            "length",
            ensure_positive=True
        )

        tools.check_pint_quantity(
            pressure,
            "pressure",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            rupture_modulus,
            "pressure",
            ensure_positive=True
        )

        safety_factor = cls.solver(
            length=length.to_base_units().magnitude,
            width=width.to_base_units().magnitude,
            thickness=thickness.to_base_units().magnitude,
            pressure=pressure.to_base_units().magnitude,
            rupture_modulus=rupture_modulus.to_base_units().magnitude
        )

        return safety_factor

    @classmethod
    def minimum_thickness(
            cls,
            length,
            width,
            safety_factor,
            pressure,
            rupture_modulus,
            unit_registry
    ):
        """
        This function calculates the thickness of a clamped rectangular window
        which gives the desired safety factor.

        Parameters
        ----------
        length : pint quantity with length units
            Window unsupported (viewing) length
        width : pint quantity with length units
            Window unsupported (viewing) width
        safety_factor : float
            Safety factor
        pressure : pint quantity with pressure units
            Design pressure differential across window at which factor of
            safety is to be calculated
        rupture_modulus : pint quantity with pressure units
            Rupture modulus of desired window material.
        unit_registry : pint unit registry
            Keeps output consistent with parent registry, avoiding conflicts

        Returns
        -------
        thickness : pint quantity
            Window thickness
        """
        quant = unit_registry.Quantity

        tools.check_pint_quantity(
            length,
            "length",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            width,
            "length",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            pressure,
            "pressure",
            ensure_positive=True
        )
        tools.check_pint_quantity(
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

        thickness = cls.solver(
            length=length.to_base_units().magnitude,
            width=width.to_base_units().magnitude,
            safety_factor=safety_factor,
            pressure=pressure.to_base_units().magnitude,
            rupture_modulus=rupture_modulus.to_base_units().magnitude
        )

        return quant(
            thickness,
            width.to_base_units().units).to(width.units.format_babel())

    @staticmethod
    def solver(
            **kwargs
    ):
        """
        This function uses sympy to solve for a missing window measurement.
        Inputs are five keyword arguments, with the following possible values:
            length
            width
            thickness
            pressure
            rupture_modulus
            safety_factor
        All of these arguments should be floats, and dimensions should be
        consistent (handling should be done in other functions, such as
        calculate_window_sf().

        Equation from:
        https://www.crystran.co.uk/userfiles/files/
        design-of-pressure-windows.pdf

        Parameters
        ----------
        kwargs

        Returns
        -------
        missing value as a float, or NaN if the result is imaginary
        """

        # Ensure that 5 keyword arguments were given
        if kwargs.__len__() != 5:
            raise ValueError("\nIncorrect number of arguments sent to solver")

        # Ensure all keyword arguments are correct
        good_arguments = [
            "length",
            "width",
            "thickness",
            "pressure",
            "rupture_modulus",
            "safety_factor"
        ]
        bad_args = []
        for arg in kwargs:
            if arg not in good_arguments:
                bad_args.append(arg)

        if len(bad_args) > 0:
            error_string = "\nBad keyword argument:"
            for arg in bad_args:
                error_string += "\n" + arg

            raise ValueError(error_string)

        # Define equation to be solved
        k_factor = 0.75  # clamped window factor
        argument_symbols = {
            "length": "var_l",
            "width": "var_w",
            "thickness": "var_t",
            "pressure": "var_p",
            "rupture_modulus": "var_m",
            "safety_factor": "var_sf"
        }
        var_l = sp.Symbol("var_l")
        var_w = sp.Symbol("var_w")
        var_t = sp.Symbol("var_t")
        var_p = sp.Symbol("var_p")
        var_m = sp.Symbol("var_m")
        var_sf = sp.Symbol("var_sf")
        expr = (
                var_l *
                var_w *
                sp.sqrt(
                    (
                            var_p *
                            k_factor *
                            var_sf /
                            (
                                    2 *
                                    var_m *
                                    (
                                            var_l ** 2 +
                                            var_w ** 2
                                    )
                            )
                    )
                ) - var_t
        )

        # Solve equation
        for arg in kwargs:
            expr = expr.subs(argument_symbols[arg], kwargs[arg])

        solution = sp.solve(expr)[0]

        if solution.is_real:
            return float(solution)
        else:
            warnings.warn("Window inputs resulted in imaginary solution.")
            return np.NaN

    @staticmethod
    def calculate_bolt_sfs(
            max_pressure,
            window_area,
            num_bolts,
            thread_size,
            thread_class,
            bolt_max_tensile,
            plate_max_tensile,
            engagement_length,
            unit_registry
    ):
        """
        Calculates bolt and plate safety factors for viewing window bolts

        Parameters
        ----------
        max_pressure : pint quantity
            Pint quantity of tube maximum pressure (absolute)
        window_area : pint quantity
            Pint quantity of window area exposed to high pressure environment
        num_bolts : int
            Number of bolts used to secure each viewing window
        thread_size : str
            Size of threads to be evaluated, e.g. '1/4-20' or '1 1/2-6'
        thread_class : str
            Class of threads to be evaluated, '2' or '3'. 'A' or 'B' are
            automatically appended for internal/external threads
        bolt_max_tensile : pint quantity
            Pint quantity of bolt (ext. thread) tensile failure stress
        plate_max_tensile : pint quantity
            Pint quantity of plate (int. thread) tensile failure stress
        engagement_length : pint quantity
            Pint quantity of total thread engagement length
        unit_registry : pint unit registry
            Keeps output consistent with parent registry, avoiding conflicts

        Returns
        -------
        safety_factors : dict
            Dictionary with keys of 'bolt' and 'plate', giving factors of safety
            for window bolts and the plate that they are screwed into.
        """
        quant = unit_registry.Quantity

        tools.check_pint_quantity(
            max_pressure,
            "pressure",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            window_area,
            "area",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            bolt_max_tensile,
            "pressure",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            plate_max_tensile,
            "pressure",
            ensure_positive=True
        )
        tools.check_pint_quantity(
            engagement_length,
            "length",
            ensure_positive=True
        )

        # convert all quantities to local unit registry
        max_pressure = quant(
            max_pressure.magnitude,
            max_pressure.units.format_babel()
        )
        window_area = quant(
            window_area.magnitude,
            window_area.units.format_babel()
        )
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

        # get total force per bolt
        window_force = (
                (max_pressure - quant(1, "atm")) * window_area / num_bolts
        )

        # get stress areas
        thread = Bolt.calculate_stress_areas(
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
        ).to_base_units()
        safety_factors["plate"] = (
                plate_max_tensile / (window_force / plate_area)
        ).to_base_units()
        return safety_factors


class Tube:
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

        Returns
        -------
        max_stress : pint quantity
            Pint quantity of maximum allowable tube stress
        """

        # look up stress-temperature limits and units
        stress_limits = cls.get_pipe_stress_limits(
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
    def get_pipe_stress_limits(
            material,
            welded,
            unit_registry
    ):
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
        This function calculates the dynamic load factor by which a detonation
        tube's static analysis should be scaled in order to account for the
        tube's response to pressure transients. DLF is based on the work of
        Shepherd [1]. Since the limits of "approximately equal to" are not
        define we assume a default value of plus or minus ten percent, thus
        plus_or_minus=0.1.

        [1] Shepherd, J. E. (2009). Structural Response of Piping to
        Internal Gas Detonation. Journal of Pressure Vessel Technology,
        131(3), 031204. https://doi.org/10.1115/1.3089497

        Parameters
        ----------
        tube_id : pint.Quantity
            Inner diameter of tube
        tube_od : pint.Quantity
            Outer diameter of tube
        cj_velocity : pint.Quantity
            CJ velocity
        elastic_modulus : pint.Quantity
            Elastic modulus of tube material
        density : pint.Quantity
            Density of tube material
        poisson_ratio : float
            Poisson ratio of tube material
        plus_or_minus : float
            Defines the band about the critical velocity which is considered
            "approximately equal to" -- the default value of 0.1 means plus
            or minus ten percent.

        Returns
        -------
        dynamic_load_factor : float
            Factor by which the tube's static maximum pressure should be
            de-rated to account for transient response to detonation waves.
        """
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
            crit_velocity * (1. + plus_or_minus),
            crit_velocity * (1. - plus_or_minus)
        ]

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
            safety_factor
    ):
        """
        Calculates the maximum allowable pressure from the tube dimensions
        and stress limits using the basic longitudinal joint formula
        on page 14 of Megyesy's Pressure Vessel Handbook, 8th ed.

        Parameters
        ----------
        tube_id : pint.Quantity
            Inner diameter of tube
        tube_od : pint.Quantity
            Outer diameter of tube
        max_stress : pint.Quantity
            Maximum allowable stress in tube material (e.g. from ASME B31.1)
        safety_factor : float
            Desired safety factor

        Returns
        -------
        max_pressure : pint.Quantity
            Pressure resulting in maximum allowable stress
        """
        mean_diameter = (tube_od + tube_id) / 2.
        wall_thickness = (tube_od - tube_id) / 2.
        max_pressure = max_stress * (2 * wall_thickness) / \
            (mean_diameter * safety_factor)

        return max_pressure

    @classmethod
    def calculate_initial_pressure(
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
            unit_registry=pint.UnitRegistry(),
            error_tol=1e-4,
            max_iterations=500,
            use_multiprocessing=False
    ):
        """
        Parameters
        ----------
        tube_id : pint.Quantity
            Inner diameter of tube
        tube_od : pint.Quantity
            Outer diameter of tube
            Internal diameter of the detonation tube
        initial_temperature : pint quantity
            Mixture initial temperature
        species_dict : dict
            Dictionary containing the species in the mixture as keys, with total
            moles or mole fractions as values
        mechanism : str
            Mechanism file name for Cantera
        max_pressure : pint quantity
            Maximum allowable pressure within the tube
        elastic_modulus : pint.Quantity
            Elastic modulus of tube material
        density : pint.Quantity
            Density of tube material
        poisson_ratio : float
            Poisson ratio of tube material
        plus_or_minus : float
            Defines the band about the critical velocity which is considered
            "approximately equal to" -- the default value of 0.1 means plus
            or minus ten percent.
        unit_registry : pint unit registry
            Unit registry for managing units to prevent conflicts with parent
            unit registry
        error_tol : float
            Relative error tolerance below which initial pressure calculations
            are considered 'good enough'
        max_iterations : int
            Maximum number of loop iterations before exit, defaults to 500
        use_multiprocessing : bool
            True to use multiprocessing. Defaults to False.

        Returns
        -------
        initial_pressure : pint quantity
            Initial mixture pressure corresponding to the tube's maximum
            allowable pressure.
        """
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

            # calculate new error, accounting for dynamic load factor
            dlf = cls.dynamic_load_factor(
                tube_id,
                tube_od,
                state["cj"]["speed"],
                elastic_modulus,
                density,
                poisson_ratio,
                plus_or_minus
            )
            error = (state["reflected"]["state"].P * dlf -
                     max_pressure.magnitude) / max_pressure.magnitude

        return initial_pressure

    @staticmethod
    def get_available_pipe_schedules(pipe_size):
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
            pipe_schedule
    ):
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
            inner_diameter=_Q(inner_diameter, "in"),
            outer_diameter=_Q(outer_diameter, "in"),
            wall_thickness=_Q(wall_thickness, "in"),
        )


class Flange:
    available_materials = TUBE_MATERIALS.Grade[
        pd.notna(TUBE_MATERIALS.Group)
    ].values

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
        Finds the minimum allowable flange class per ASME B16.5 for a give
        flange temperature and tube pressure.

        Returns
        -------
        flange_class: str
            String representing the minimum allowable flange class
        """
        max_pressure = tools.parse_quant_input(
            max_pressure,
            unit_registry
        )
        temperature = tools.parse_quant_input(
            temperature,
            unit_registry
        )
        tools.check_pint_quantity(
            max_pressure,
            "pressure",
            True
        )
        tools.check_pint_quantity(
            temperature,
            "temperature",
            True
        )
        max_pressure = tools.parse_quant_input(
            max_pressure,
            _U
        )
        temperature = tools.parse_quant_input(
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
            raise ValueError("\nPressure out of range.")

        # ensure temperature is within bounds
        if (
                temperature < df_limits["Temperature"].min()
        ) or (
                temperature > df_limits["Temperature"].max()
        ):
            raise ValueError("\nTemperature out of range.")

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
            ser_result > max_pressure.to_base_units().magnitude
        ][0]

    @staticmethod
    def get_max_pressure(
            flange_class,
            temperature,
            material
    ):
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

        temperature = tools.parse_quant_input(temperature, _U)

        df_limits = FLANGE_LIMITS[group]
        # ensure temperature is in range
        if (
                temperature < df_limits["Temperature"].min()
        ) or (
                temperature > df_limits["Temperature"].max()
        ):
            raise ValueError("\nTemperature out of range.")

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

        pressure = tools.parse_quant_input(pressure, _U)

        df_limits = FLANGE_LIMITS[group]
        # ensure pressure is in range
        if (
                pressure < df_limits[flange_class].min()
        ) or (
                pressure > df_limits[flange_class].max()
        ):
            raise ValueError("\nPressure out of range.")

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