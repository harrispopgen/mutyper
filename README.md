mutyper
====

A Python package for annotating the local ancestral sequence context of biallelic SNPs, featuring

- CLI for integration in pipelines analyzing VCF/BCF data
- Python API

Installation
---
1. clone the repository
```bash
$ git clone https://github.com/harrispopgen/mutyper
```
2. install with pip
```bash
$ pip install mutyper
```
 **Developer note:** for editable installation, instead use
```bash
$ pip install -e mutyper
```
3. To run the demonstration Jupyter notebook [`demo.ipynb`](demo.ipynb), you'll also need to install Jupyter, bcftools, and tabix. If you're working in a Conda environment, this can be done with
```bash
$ conda install bcftools tabix jupyter --channel conda-forge --channel bioconda
```

Command line usage
---
- ### main help and subcommand list
```bash
$ mutyper -h
```

- ### help for each subcommand
```bash
$ mutyper <subcommand> -h
```

Python API
---

- ### `ancestor` module
```python
>>> from mutyper import ancestor
>>> help(ancestor)
```
