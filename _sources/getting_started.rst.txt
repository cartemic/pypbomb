Installation
------------

pybomb requires a few dependencies in order to operate. They are

* `Cantera <https://cantera.org/>`_

* `NumPy <https://numpy.org/>`_

* `Pandas <https://pandas.pydata.org/>`_

* `Pint <https://pint.readthedocs.io/en/stable/>`_, and

* `SciPy <https://scipy.org/>`_.

Note that at least a couple of these can be installed via multiple sources
(e.g. pip, conda, your favorite package manager). If you install them manually
your desired python interpreter must be able to access them, and your best bet
for smooth operation is for them to be the versions specified in
requirements.txt. The easiest option is to install pypbomb as described below
and let pip sort it out for you.

In addition the the listed dependencies, the ``sd`` portion of this package
contains a small portion of
`SDToolbox <https://shepherd.caltech.edu/EDL/PublicResources/sdt/>`_ adapted
for use here in accordance with the licensing agreement and with permission
from the authors. If you are planning to do more general detonation simulation
work, it is a great resource.

For Use
^^^^^^^

In order to use pypbomb, perform the following steps.

#. Download and extracted the latest release from
   `github <https://github.com/cartemic/pypbomb/releases/latest>`_ and navigate
   to the top level directory of pypbomb (e.g. ``~/Downloads/pypbomb/``) in your
   terminal of choice.

#. Install pypbomb using pip.

   .. code-block:: bash

      $ pip install .

#. Delete the pypbomb directory that you downloaded, if desired.

For Development
^^^^^^^^^^^^^^^

If you would like to modify pypbomb, the recommended method is to install it for
development in a conda environment. For direction on creating a conda
environment, see the
`conda documentation <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.

#. Move the downloade pypbomb directory to where you want it to live
   permanently.

   .. code-block:: bash

      $ mkdir ~/my/favorite/dir/
      $ mv pypbomb/ ~/my/favorite/dir/

#. Navigate to that location in the terminal of your choice.

   .. code-block:: bash

      $ cd ~/my/favorite/dir/pypbomb

#. Activate your conda environment and install the package requirements.

   .. code-block:: bash

      $ source activate myenv
      (myenv)$ pip install -r requirements.txt

#. Install pypbomb in development mode.

   .. code-block:: bash

      (myenv)$ conda develop .

In this case, *do not* delete the pypbomb directory, since this is where pypbomb
lives now. The source code is in pypbomb/pypbomb/.

If you are planning to make changes to pypbomb, you may also wish to make
changes to the documentation. The documentation is built using sphinx,
and uses the readthedocs theme. Its source is located within pypbomb/docs/.
If you've just installed pypbomb, now would be a good time to install the
documentation dependencies.

.. code-block:: bash

   (myenv)$ pip install -r docs/requirements.txt

Once that is out of the way, you are free to modify the files at will. To
build them locally, run the following command from within the docs/ directory.

.. code-block:: bash

   (myenv)$ make clean html

Sweet! Your new documentation should now be in pypbomb/docs/_build/html.

If you've made and documented some modifications that you'd like to contribute,
please see the `section on contributing <#contributing>`_.

Testing
-------

To test this code, you must first install the testing requirement (if you don't
have them already).

.. code-block:: bash

   (myenv)$ pip install -r requirements_test.txt

Once you have the requirements installed you can test pypbomb and check the
coverage from the top level pypbomb directory.

.. code-block:: bash

   (myenv)$ pytest -vv cov=./

Contributing
------------
Pypbomb is open source, and contributing is encouraged. If you would like to
contribute some modifications you've made, please submit a pull request to
the development branch.

This project uses isort, black, and flake8 to maintain some semblance of code
quality. It is my experience that this is helpful for further development and
maintenance. Before submitting your code for PR, please use these tools first
in order to maintain a somewhat consistent code style. This can be easily
done using the ``prepare.sh`` script in the top level project directory.

In addition to running the script manually, it is recommended to install a
pre-commit hook in order to prevent any commits that aren't properly styled or
have linter gripes. This just makes it easier on everybody. pre-commit is one
of the dev requirements, so if you've installed those, all you need to do is
run ``pre-commit install`` and you should be all set up!
