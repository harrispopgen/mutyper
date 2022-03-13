Installation
############

Pip
===

.. code-block:: bash

  pip install mutyper


Conda
=====

.. todo::
  setup Conda install

Additional dependencies
=======================

If installation fails with a build failure from ``cyvcf2``, you may need to install a few C libraries.
To create and activate a Conda environment with these, e.g. named ``mutyperenv``:

.. code-block:: bash

  conda create -n mutyperenv -c anaconda libcurl bzip2 pip
  conda activate mutyperenv

Then ``mutyper`` can be installed using pip as shown above.
If you do not already have a C compiler (``gcc``), you can append one in the ``conda create ...`` command above (``gcc_linux-64`` for linux, ``clang_osx-64`` for macos).
