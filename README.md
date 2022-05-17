# pypbomb

[![Black badge](https://camo.githubusercontent.com/d91ed7ac7abbd5a6102cbe988dd8e9ac21bde0a73d97be7603b891ad08ce3479/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f636f64652532307374796c652d626c61636b2d3030303030302e737667)](https://github.com/psf/black)
[![Isort badge](https://camo.githubusercontent.com/fe4a658dd745f746410f961ae45d44355db1cc0e4c09c7877d265c1380248943/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f253230696d706f7274732d69736f72742d2532333136373462313f7374796c653d666c6174266c6162656c436f6c6f723d656638333336)](https://github.com/PyCQA/isort)
[![Travis build badge](https://api.travis-ci.com/cartemic/pypbomb.svg?branch=master)](https://travis-ci.com/github/cartemic/pypbomb)
[![Codecov badge](https://codecov.io/gh/cartemic/pypbomb/branch/master/graph/badge.svg)](https://codecov.io/gh/cartemic/pypbomb)
[![Documentation badge](https://github.com/cartemic/pypbomb/workflows/Documentation/badge.svg)](https://cartemic.github.io/pypbomb/)
[![status](https://joss.theoj.org/papers/288143539a2e4f5ab1b83cbee4e17052/status.svg)](https://joss.theoj.org/papers/288143539a2e4f5ab1b83cbee4e17052)

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



