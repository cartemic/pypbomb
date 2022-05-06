Installation and Testing
========================

Installation
------------

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
please feel free to submit a pull request to the develop branch of this repo.

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
