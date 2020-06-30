Ancestral :math:`k`-mer Mutation Types for SNP Data
###################################################

``mutyper`` is a Python package for augmenting genomic variation data in VCF/BCF format with information about the local ancestral sequence context of each variant.
Studying mutation spectra usually begins with computing the mutation type for SNPs in a VCF/BCF file
To polarize SNPs from REF/ALT to ancestral/derived and determine their local context, we need an estimate of the ancestral sequence, which usually takes the form of a FASTA file.
``mutyper`` is a Python package for doing this.
It also has a command-line interface for easily integrating in more complex bioinformatic pipelines.

.. toctree::
  :maxdepth: 1
  :caption: User Guide

  install
  notebooks/quickstart
  cite

.. toctree::
   :maxdepth: 3
   :caption: CLI Documentation

   cmd

.. toctree::
   :maxdepth: 3
   :caption: API Documentation

   mutyper

.. toctree::
  :maxdepth: 1
  :caption: Developer Documentation

  developer

.. toctree::
   :maxdepth: 1
   :caption: Notes

   CHANGELOG
   faq


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
