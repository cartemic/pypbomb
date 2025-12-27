import sqlite3
import typing
from dataclasses import dataclass
from functools import cached_property
from typing import Optional

import numpy as np

from pypbomb._types import FloatArray
from pypbomb.data import _db, _units


def available_materials(flange_class: Optional[int] = None) -> set[str]:
    if flange_class is None:
        return _all_available_materials()
    else:
        return _available_materials_for_class(flange_class)


def _all_available_materials() -> set[str]:
    with _db.connection() as conn:
        result = conn.execute(
            """
                select mat.grade
                from material mat
                where mat.material_group in (
                    select distinct fla.material_group
                    from b16_flange_rating fla
                );
            """
        ).fetchall()
    return set(item[0] for item in result)


def _available_materials_for_class(flange_class: int) -> set[str]:
    with _db.connection() as conn:
        result = conn.execute(
            """
                select mat.grade
                from material mat
                where mat.material_group in (
                    select fla.material_group
                    from b16_flange_rating fla
                    where fla.flange_class = :flange_class
                );
            """,
            {"flange_class": flange_class},
        ).fetchall()
    return set(item[0] for item in result)


def available_classes(material: str) -> set[int]:
    with _db.connection() as conn:
        result = conn.execute(
            """
            select distinct fla.flange_class
            from b16_flange_rating fla
            where fla.material_group in (
                select mat.material_group
                from material mat
                where mat.grade = :material
            );
            """,
            {"material": material},
        ).fetchall()
    return set(item[0] for item in result)


@dataclass(frozen=True)
class FlangePressureLimits:
    temperature: FloatArray
    pressure_limit: FloatArray

    @cached_property
    def _temp_min(self):
        return self.temperature.min()

    @cached_property
    def _temp_max(self):
        return self.temperature.max()

    @cached_property
    def _press_min(self):
        return self.pressure_limit.min()

    @cached_property
    def _press_max(self):
        return self.pressure_limit.max()

    def max_pressure(self, temperature: float) -> float:
        if not self._temp_min <= temperature <= self._temp_max:
            raise ValueError(f"Pressure must be {self._temp_min} <= temperature <= {self._temp_max}")
        return typing.cast(float, np.interp(x=temperature, xp=self.temperature, fp=self.pressure_limit))

    def max_temperature(self, pressure: float) -> float:
        if not self._press_min <= pressure <= self._press_max:
            raise ValueError(f"Pressure must be {self._press_min} <= pressure <= {self._press_max}")
        # Limits are sorted by ascending temperature, meaning pressure limit is monotonically decreasing.
        # Therefore, we need to flip the arrays before interpolating.
        return typing.cast(
            float,
            np.interp(x=pressure, xp=self.pressure_limit[::-1], fp=self.temperature[::-1]),
        )


def pressure_limits_opt(material: str, flange_class: int) -> Optional[FlangePressureLimits]:
    with _db.connection() as conn:
        conn.row_factory = sqlite3.Row
        results = conn.execute(
            """
            select
                fla.temperature_f,
                fla.max_pressure_kpsi
            from b16_flange_rating fla
            where fla.material_group in (
                select mat.material_group
                from material mat
                where
                    mat.grade = :material
                    and fla.flange_class = :flange_class
            )
            order by fla.temperature_f asc;
            """,
            {"material": material, "flange_class": flange_class},
        ).fetchall()
        n_entries = len(results)
        if n_entries == 0:
            return None
        temperature = np.zeros(n_entries) * np.nan
        pressure_limit = np.zeros(n_entries) * np.nan
        for i, entry in enumerate(results):
            temperature[i] = entry["temperature_f"]
            pressure_limit[i] = entry["max_pressure_kpsi"]
        return FlangePressureLimits(
            temperature=_units.f_to_k(temperature),
            pressure_limit=_units.kpsi_to_pa(pressure_limit),
        )


def pressure_limits(material: str, flange_class: int) -> FlangePressureLimits:
    maybe_pressure = pressure_limits_opt(material, flange_class)
    if maybe_pressure is None:
        raise ValueError(f"Unable to find pressure limits for class {flange_class} flanges made of {material}")
    return maybe_pressure


def minimum_flange_class_opt(temperature: float, pressure: float, material: str) -> Optional[int]:
    with _db.connection() as conn:
        flange_classes = conn.execute(
            "select distinct flange_class from b16_flange_rating order by flange_class asc"
        ).fetchall()
        conn.row_factory = sqlite3.Row
        for (flange_class,) in flange_classes:
            results = conn.execute(
                """
                select
                    fla.flange_class,
                    fla.temperature_f,
                    fla.max_pressure_kpsi
                from b16_flange_rating fla
                where fla.material_group in (
                    select mat.material_group
                    from material mat
                    where
                        mat.grade = :material
                        and fla.flange_class = :flange_class
                )
                """,
                {"material": material, "flange_class": flange_class},
            ).fetchall()
            n_entries = len(results)
            temp_limit = np.zeros(n_entries) * np.nan
            pressure_limit = np.zeros(n_entries) * np.nan
            for i, entry in enumerate(results):
                temp_limit[i] = entry["temperature_f"]
                pressure_limit[i] = entry["max_pressure_kpsi"]

            temp_limit = _units.f_to_k(temp_limit)
            pressure_limit = _units.kpsi_to_pa(pressure_limit)

            if not temp_limit.min() <= temperature <= temp_limit.max():
                continue
            elif pressure > pressure_limit.max():
                continue
            elif pressure <= np.interp(x=temperature, xp=temp_limit, fp=pressure_limit):
                return flange_class

    return None


def minimum_flange_class(temperature: float, pressure: float, material: str) -> float:
    maybe_class = minimum_flange_class_opt(temperature, pressure, material)
    if maybe_class is None:
        raise ValueError(
            f"Unable to find suitable flange class for {material} at P = {pressure} Pa and T = {temperature} K"
        )
    return maybe_class
