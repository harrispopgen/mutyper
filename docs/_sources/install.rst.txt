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

If installation fails with a build failure from ``cyvcf2``, you may need to install libcurl:

**linux:**

.. code-block:: bash

  sudo apt-get update; sudo apt-get install -y libcurl4-openssl-dev
