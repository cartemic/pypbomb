def temperature(t: float):
    if t < 0:
        raise ValueError(f"Nonphysical temperature: {t}. Must be > 0.")


def pressure(p: float):
    if p < 0:
        raise ValueError(f"Nonphysical pressure: {p}. Must be > 0.")


def length(l: float):
    if l < 0:
        raise ValueError(f"Nonphysical length: {l}. Must be > 0.")


def area(p: float):
    if p < 0:
        raise ValueError(f"Nonphysical area: {p}. Must be > 0.")


def rupture_modulus(r: float):
    if r < 0:
        raise ValueError(f"Nonphysical rupture modulus: {r}. Must be > 0.")
