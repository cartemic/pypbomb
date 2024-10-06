import warnings
from dataclasses import asdict, dataclass
from typing import Optional

import numpy as np

from . import _db, _units


def available_sizes() -> set[str]:
    """
    :return: Available nominal pipe sizes
    """
    with _db.connection() as conn:
        return set(
            result[0]
            for result in conn.execute(
                """
            select pipe_size
            from pipe_od
            where pipe_size in (
                select pipe_size from pipe_thk
            );
            """
            ).fetchall()
        )


def available_schedules(size: str) -> set[str]:
    """
    :param size: Nominal pipe size, e.g. "2 1/2"
    :return: Available pipe schedules for ``size``
    """
    with _db.connection() as conn:
        return set(
            result[0]
            for result in conn.execute(
                """
            select distinct schedule
            from pipe_thk
            where pipe_size=:size;
            """,
                {"size": size},
            ).fetchall()
        )


def outer_diameter_opt(size: str) -> Optional[float]:
    """
    :param size: Nominal pipe size, e.g. "2 1/2"
    :return: Pipe outer diameter (m)
    """
    with _db.connection() as conn:
        result = conn.execute(
            """
            select od_in
            from pipe_od
            where
                pipe_size=:size
                AND od_in > 0;
            """,
            {"size": size},
        ).fetchone()
    if result is None:
        return None
    return _units.in_to_m(result[0])


def outer_diameter(size: str) -> float:
    """
    :param size: Nominal pipe size, e.g. "2 1/2"
    :return: Pipe outer diameter (m)
    """
    maybe_diameter = outer_diameter_opt(size)
    if maybe_diameter is None:
        raise ValueError(f"No outer diameter found for pipe size {size}")
    return maybe_diameter


def wall_thickness_opt(size: str, schedule: str) -> Optional[float]:
    """
    :param size: Nominal pipe size, e.g. "2 1/2"
    :param schedule: Pipe schedule, e.g. "XXH"
    :return: Pipe wall thickness (m)
    """
    with _db.connection() as conn:
        result = conn.execute(
            """
            select thk_in
            from pipe_thk
            where
                pipe_size=:size
                AND schedule=:schedule
                AND thk_in > 0;
            """,
            {
                "size": size,
                "schedule": schedule,
            },
        ).fetchone()
    if result is None:
        return None
    return _units.in_to_m(result[0])


def wall_thickness(size: str, schedule: str) -> float:
    """
    :param size: Nominal pipe size, e.g. "2 1/2"
    :param schedule: Pipe schedule, e.g. "XXH"
    :return: Pipe wall thickness (m)
    """
    maybe_thickness = wall_thickness_opt(size, schedule)
    if maybe_thickness is None:
        raise ValueError(f"No wall thickness found for pipe size {size} with schedule {schedule}")
    return maybe_thickness


def inner_diameter_opt(outer: Optional[float], thk: Optional[float]) -> Optional[float]:
    if None in (outer, thk):
        return None
    elif thk >= outer / 2:
        warnings.warn("Insufficient wall thickness for given diameter")
        return None
    return outer - 2 * thk


def inner_diameter(outer: float, thk: float) -> float:
    maybe_diameter = inner_diameter_opt(outer, thk)
    if outer_diameter is None:
        raise ValueError("Invalid pipe dimensions encountered, please fix lookup table")
    return maybe_diameter


def mean_diameter_opt(outer: Optional[float], inner: Optional[float]) -> Optional[float]:
    if None in (outer, inner):
        return None
    elif inner >= outer:
        warnings.warn("Outer diameter must be greater than inner diameter")
        return None
    return (outer + inner) / 2


def mean_diameter(outer: float, inner: float) -> float:
    maybe_diameter = mean_diameter_opt(outer, inner)
    if outer_diameter is None:
        raise ValueError("Invalid pipe dimensions encountered, please fix lookup table")
    return maybe_diameter


@dataclass(frozen=True)
class PipeDimensions:
    """
    Pipe dimensions (m)
    """

    inner_diameter: float
    "Inner Diameter (m)"
    outer_diameter: float
    "Outer Diameter (m)"
    mean_diameter: float
    "Mean Diameter (m)"
    wall_thickness: float
    "Wall Thickness (m)"

    # Not sure why PyCharm is griping about PipeDimensions not being DataclassInstance
    # noinspection PyTypeChecker
    def __eq__(self, other: "PipeDimensions") -> bool:
        other_dict = asdict(other)
        return all(np.isclose(v, other_dict[k]) for k, v in asdict(self).items())


def dimensions_opt(size: str, schedule: str) -> Optional[PipeDimensions]:
    maybe_outer = outer_diameter_opt(size)
    maybe_thk = wall_thickness_opt(size, schedule)
    maybe_inner = inner_diameter_opt(maybe_outer, maybe_thk)
    maybe_mean = mean_diameter_opt(maybe_outer, maybe_inner)
    if None in (maybe_outer, maybe_thk, maybe_inner):
        return None
    return PipeDimensions(
        inner_diameter=maybe_inner,
        outer_diameter=maybe_outer,
        mean_diameter=maybe_mean,
        wall_thickness=maybe_thk,
    )


def dimensions(size: str, schedule: str) -> PipeDimensions:
    outer = outer_diameter(size)
    thk = wall_thickness(size, schedule)
    inner = inner_diameter(outer, thk)
    mean = mean_diameter(outer, inner)
    return PipeDimensions(
        inner_diameter=inner,
        outer_diameter=outer,
        mean_diameter=mean,
        wall_thickness=thk,
    )
