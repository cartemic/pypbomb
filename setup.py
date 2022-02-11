import sys
from os import path, system

from setuptools import setup

here = path.abspath(path.dirname(__file__))

# conda build is being a real jerk and claiming I don't have six (I do), so
# we're going to do some hack-y stuff here. Force pint and cantera to conda
# install and just... roll with it.
if system("conda list"):
    # nonzero return means failure... pls do the conda
    # IDK how to test this yet so I'm just going to hope it works until told
    # otherwise.
    system("pip install conda")

try:
    import cantera
except ModuleNotFoundError:
    system("conda install -c cantera cantera -y")

good_pint_version = "0.11"
try:
    import pint
    if pint.__version__ != good_pint_version:
        raise ModuleNotFoundError
except ModuleNotFoundError:
    system(f"conda install -c conda-forge pint={good_pint_version} -y")

with open(path.join(here, "pypbomb", "_version.py")) as f:
    __version__ = ""
    exec(f.read())

with open(path.join(here, "README.md")) as f:
    readme = f.read()

with open(path.join(here, "CHANGELOG.md")) as f:
    changelog = f.read()


desc = readme + "\n\n" + changelog
try:
    import pypandoc
    long_description = pypandoc.convert_text(desc, 'rst', format='md')
    with open(path.join(here, 'README.rst'), 'w') as rst_readme:
        rst_readme.write(long_description)
except (ImportError, OSError, IOError):
    long_description = desc

install_requires = [
    # "cantera",  # we're not doing conda install anymore. is broken.
    "numpy",
    "pytest",
    "pandas",
    # "pint",     # we're not doing conda install anymore. is broken.
    "six"
]
tests_require = [
    "pytest",
    "pytest-cov",
    "mock",
]

needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
setup_requires = ["pytest-runner"] if needs_pytest else []

setup(
    name="pypbomb",
    version=__version__,
    python_requires=">=3.7.*",
    packages=["pypbomb", "pypbomb.tests"],
    url="https://github.com/cartemic/pypbomb",
    license="BSD-3",
    author="Mick Carter",
    author_email="cartemic@oregonstate.edu",
    package_dir={"pypbomb": "pypbomb"},
    package_data={
        "pypbomb": [path.join("lookup_data", "*")],
        "pypbomb.tests": [path.join("lookup_data", "*")]
    },
    description="Tools for designing a detonation tube",
    long_description=long_description,
    install_requires=install_requires,
    tests_require=tests_require,
    setup_requires=setup_requires,
    py_modules=["six"],
)
