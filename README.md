# pypbomb

[![Black badge](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Isort badge](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat)](https://github.com/PyCQA/isort)
[![Travis build badge](https://api.travis-ci.com/cartemic/pypbomb.svg?branch=master)](https://travis-ci.com/github/cartemic/pypbomb)
[![Codecov badge](https://codecov.io/gh/cartemic/pypbomb/branch/master/graph/badge.svg)](https://codecov.io/gh/cartemic/pypbomb)
[![Documentation badge](https://github.com/cartemic/pypbomb/workflows/Documentation/badge.svg)](https://cartemic.github.io/pypbomb/)
[![status](https://joss.theoj.org/papers/288143539a2e4f5ab1b83cbee4e17052/status.svg)](https://joss.theoj.org/papers/288143539a2e4f5ab1b83cbee4e17052)
[![DOI](https://zenodo.org/badge/130134402.svg)](https://zenodo.org/badge/latestdoi/130134402)

---

For when you want the detonation to stay on the inside.

## Overview

This package was made to help with the design of a closed-end detonation tube, and was written as part of my PhD research in the [propulsion lab](http://research.engr.oregonstate.edu/blunckgroup/propulsion-laboratory) at Oregon Sate University. It contains some tools that can be used to iteratively determine initial design parameters for a closed-end detonation tube.

The main items of interest are:

* `Bolt`: Methods relating to bolt calculations and property lookup.
* `DDT`: Methods for estimating the deflagration-to-detonation transition (DDT).
* `Flange`: Methods pertaining to flange class and max pressure calculations.
* `Tube`: Methods for designing and determining operational limits of a detonation tube.
* `Window`: Methods for the design of viewing windows for optical access.

Full documentation can be found [here](https://cartemic.github.io/pypbomb/).

## Installation

Once cloned to your local machine, `pypbomb` can be installed via pip from within the project's top-level directory:

```bash
pip install .
```



Further instructions can be found in the full documentation if needed.

## Short Example

To get an estimate for the maximum initial pressure that a tube with a given geometry can handle:

```python
import pint
from pypbomb import Tube

ureg = pint.UnitRegistry()
quant = ureg.Quantity

initial_pressure = Tube.calculate_max_initial_pressure(
    tube_id=quant(5.76, "in"),
    tube_od=quant(6.625, "in"),
    initial_temperature=quant(70, "degF"),
    species_dict={"H2":2, "O2":1},
    mechanism="gri30.yaml",
    max_pressure=(2, "kpsi"),
    elastic_modulus=quant(200, "GPa"),
    density=quant(8, "g/cm**3"),
    poisson_ratio=0.27,
)
```



