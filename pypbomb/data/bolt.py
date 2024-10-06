import sqlite3
from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from string import Template
from typing import Optional

from pypbomb.data import _db, _units


@dataclass(frozen=True)
class DiameterRange:
    """
    Diameter range (m)
    """

    min: float
    max: float


@dataclass(frozen=True)
class InternalDiameterAnsi:
    """
    ANSI internal diameters, converted to meters
    """

    basic: float
    minor: DiameterRange
    pitch: DiameterRange
    major: float


@dataclass(frozen=True)
class InternalThreadAnsi:
    """
    ANSI internal thread specs, converted to meters
    """

    thread_size: str
    thread_class: str
    threads_per_inch: int
    diameter: InternalDiameterAnsi

    @cached_property
    def threads_per_meter(self) -> float:
        return _units.tpi_to_tpm(self.threads_per_inch)


@dataclass(frozen=True)
class ExternalDiameterAnsi:
    """
    ANSI external diameters, converted to meters
    """

    basic: float
    major: DiameterRange
    pitch: DiameterRange
    minor: float


@dataclass(frozen=True)
class ExternalThreadAnsi:
    """
    ANSI external thread specs, converted to meters
    """

    thread_size: str
    thread_class: str
    threads_per_inch: int
    diameter: ExternalDiameterAnsi

    @cached_property
    def threads_per_meter(self) -> float:
        return _units.tpi_to_tpm(self.threads_per_inch)


class ThreadType(Enum):
    internal = "internal"
    external = "external"

    def __str__(self):
        return self.value


def available_thread_sizes(thread_type: ThreadType, thread_class: Optional[str] = None) -> set[str]:
    if thread_type not in ThreadType:
        return set()
    query = Template("select distinct thread_size from ansi_inch_thread_${thread_type} ${maybe_where};")
    maybe_where = "where thread_class=:thread_class" if thread_class is not None else ""
    with _db.connection() as conn:
        return set(
            result[0]
            for result in conn.execute(
                query.substitute({"thread_type": thread_type, "maybe_where": maybe_where}),
                {"thread_class": thread_class or ""},
            ).fetchall()
        )


def available_thread_classes(thread_type: ThreadType, thread_size: Optional[str] = None) -> set[str]:
    if thread_type not in ThreadType:
        return set()
    query = Template("select distinct thread_class from ansi_inch_thread_${thread_type} ${maybe_where};")
    maybe_where = "where thread_size=:thread_size" if thread_size is not None else ""
    with _db.connection() as conn:
        return set(
            result[0]
            for result in conn.execute(
                query.substitute({"thread_type": thread_type, "maybe_where": maybe_where}),
                {"thread_size": thread_size or ""},
            ).fetchall()
        )


def properties_opt(
    thread_type: ThreadType,
    *,
    thread_size: str,
    thread_class: str,
) -> Optional[InternalThreadAnsi | ExternalThreadAnsi]:
    if thread_type is ThreadType.internal:
        return internal_thread_properties_opt(thread_class=thread_class, thread_size=thread_size)
    elif thread_type is ThreadType.external:
        return external_thread_dimensions(thread_class, thread_size)
    else:
        return None


def internal_thread_properties_opt(*, thread_class: str, thread_size: str) -> Optional[InternalThreadAnsi]:
    with _db.connection() as conn:
        conn.row_factory = sqlite3.Row
        result = conn.execute(
            """
            select
                d_basic_in,
                d_minor_min_in,
                d_minor_max_in,
                d_pitch_min_in,
                d_pitch_max_in,
                d_major_in
            from ansi_inch_thread_internal
            where
                thread_class=:thread_class
                AND thread_size=:thread_size
            ;
            """,
            {"thread_class": thread_class, "thread_size": thread_size},
        ).fetchone()
    if result is None:
        return result
    return InternalThreadAnsi(
        thread_size=thread_size,
        thread_class=thread_class,
        threads_per_inch=_threads_per_inch(thread_size),
        diameter=InternalDiameterAnsi(
            basic=_units.in_to_m(result["d_basic_in"]),
            minor=DiameterRange(
                min=_units.in_to_m(result["d_minor_min_in"]),
                max=_units.in_to_m(result["d_minor_max_in"]),
            ),
            pitch=DiameterRange(
                min=_units.in_to_m(result["d_pitch_min_in"]),
                max=_units.in_to_m(result["d_pitch_max_in"]),
            ),
            major=_units.in_to_m(result["d_major_in"]),
        ),
    )


def external_thread_properties_opt(*, thread_class: str, thread_size: str) -> Optional[ExternalThreadAnsi]:
    with _db.connection() as conn:
        conn.row_factory = sqlite3.Row
        result = conn.execute(
            """
            select
                d_basic_in,
                d_major_min_in,
                d_major_max_in,
                d_pitch_min_in,
                d_pitch_max_in,
                d_minor_in
            from ansi_inch_thread_external
            where
                thread_class=:thread_class
                AND thread_size=:thread_size
            ;
            """,
            {"thread_class": thread_class, "thread_size": thread_size},
        ).fetchone()
    if result is None:
        return result
    return ExternalThreadAnsi(
        thread_size=thread_size,
        thread_class=thread_class,
        threads_per_inch=_threads_per_inch(thread_size),
        diameter=ExternalDiameterAnsi(
            basic=_units.in_to_m(result["d_basic_in"]),
            major=DiameterRange(
                min=_units.in_to_m(result["d_major_min_in"]),
                max=_units.in_to_m(result["d_major_max_in"]),
            ),
            pitch=DiameterRange(
                min=_units.in_to_m(result["d_pitch_min_in"]),
                max=_units.in_to_m(result["d_pitch_max_in"]),
            ),
            minor=_units.in_to_m(result["d_minor_in"]),
        ),
    )


def _threads_per_inch(thread_size: str) -> int:
    try:
        return int(thread_size.split("-")[-1])
    except ValueError:
        raise ValueError(f"Unable to parse TPI from thread size: {thread_size}")


def dimensions(
    thread_type: ThreadType,
    *,
    thread_size: str,
    thread_class: str,
) -> InternalThreadAnsi | ExternalThreadAnsi:
    maybe_properties = properties_opt(thread_type, thread_class=thread_class, thread_size=thread_size)
    if maybe_properties is None:
        raise ValueError(f"No {thread_type} thread properties found for class {thread_class}, size {thread_size}")
    return maybe_properties


def internal_thread_dimensions(thread_class: str, thread_size: str) -> InternalThreadAnsi:
    maybe_properties = internal_thread_properties_opt(thread_class=thread_class, thread_size=thread_size)
    if maybe_properties is None:
        raise ValueError(f"No internal thread properties found for class {thread_class}, size {thread_size}")
    return maybe_properties


def external_thread_dimensions(thread_class: str, thread_size: str) -> ExternalThreadAnsi:
    maybe_properties = external_thread_properties_opt(thread_class=thread_class, thread_size=thread_size)
    if maybe_properties is None:
        raise ValueError(f"No external thread properties found for class {thread_class}, size {thread_size}")
    return maybe_properties
