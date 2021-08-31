Open source code repository
===========================

All code is freely available at `<https://github.com/harrispopgen/mutyper>`_

Developer tools
===============

Install and activate Conda developer environment::

  conda env create -f env.yml
  conda activate mutyperdocs

Local install of ``mutyper``::

  pip install -e .

Run tests::

  make test

Format code::

  make format

Lint::

  make lint

Build docs locally (you can then see the generated documentation in ``docsrc/_build/html/index.html``.)::

  make docs

Deploy docs to ``docs/`` directory for GitHub Pages::

  make deploy


Todo list
=========

.. todolist::
