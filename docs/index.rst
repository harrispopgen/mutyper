Ancestral :math:`k`-mer Mutation Types for SNP Data
###################################################

``mutyper`` is a command line utility and Python package for augmenting genomic variation data in VCF/BCF format with information about the local ancestral sequence context of each variant.
Studying mutation spectra usually begins with computing the mutation type for SNPs in a VCF/BCF file.
To polarize SNPs from REF/ALT to ancestral/derived and determine their local context, we need an estimate of the ancestral sequence, which usually takes the form of a FASTA file.
The ``mutyper`` command-line interface makes it easy to include such analysis in more complex bioinformatic pipelines.
The Python API can be used for exploratory mutation type analysis in Python scripts or Jupyter notebooks.

.. toctree::
  :maxdepth: 1
  :caption: User Guide

  install
  quickstart
  cite

.. toctree::
   :maxdepth: 3
   :caption: CLI

   cmd

.. toctree::
   :maxdepth: 3
   :caption: API

   mutyper

.. toctree::
  :maxdepth: 1
  :caption: Developer

  developer

.. toctree::
   :maxdepth: 1
   :caption: Notes

   CHANGELOG


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
