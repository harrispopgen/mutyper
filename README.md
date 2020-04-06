![](logo.png)

A Python package for annotating the local ancestral sequence context of biallelic SNPs, featuring

- CLI for integration in pipelines analyzing VCF/BCF data
- Python API

Installation
---

- Basic install with pip
```bash
$ pip install git+https://github.com/harrispopgen/mutyper
```

- Developer installation: add `-e` for editable (clones repo to `./src/mutyper`)
```bash
$ pip install -e git+https://github.com/harrispopgen/mutyper#egg=mutyper
```

- To run the demonstration Jupyter notebook [`demo.ipynb`](demo.ipynb), you'll also need to install Jupyter, bcftools, and tabix. If you're working in a Conda environment, this can be done with
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
