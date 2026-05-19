Installation of astroMUGS
************


How to obtain radprocess
========================

Clone the GitHub repository, create a virtual environment, and install in
editable mode::

    git clone https://github.com/sachagavino/astroMUGS.git
    cd astromugs
    python -m venv .venv --prompt astromugs
    source .venv/bin/activate
    python -m pip install -e .

This creates a virtual environment (``.venv/``) inside the ``astromugs/``
directory, activates it, and installs radprocess with all its Python
dependencies (numpy, scipy, zarr, matplotlib, Xarray, etc.).

.. note::

   Always use ``python -m pip`` instead of bare ``pip``. On some systems,
   ``pip`` may point to a different Python environment than the one activated
   by the venv, which leads to packages being installed in the wrong place.


Every time you open a new terminal and want to use astromugs, you need to
activate the environment first::

    cd /path/to/radprocess
    source .venv/bin/activate

To update your local copy later::

    git pull
    python -m pip install -e .
    
 


Running the code
=================

In a python script or notebook, you can import the astromugs pipeline::

    import astromugs.pipeline as pipeline



