---
title: '`mutyper`: assigning and summarizing mutation types for analyzing germline mutation spectra'
tags:
  - genomics
  - computational biology
  - bioinformatics
  - mutation spectrum
  - python
authors:
  - name: William S. DeWitt
    orcid: 0000-0002-6802-9139
    affiliation: 1
    corresponding: true
  - name: Luke Zhu
    orcid: 0000-0002-6324-1464
    affiliation: 2
  - name: Mitchell R. Vollger
    orcid: 0000-0002-8651-1615
    affiliation: 2
  - name: Michael E. Goldberg
    orcid: 0000-0003-3310-6349
    affiliation: 3
  - name: Andrea Talenti
    orcid: 0000-0003-1309-3667
    affiliation: 4
  - name: Annabel Beichman
    orcid: 0000-0002-6991-587X
    affiliation: 2
  - name: Kelley Harris
    orcid: 0000-0003-0302-2523
    affiliation: 2
affiliations:
 - name: Department of Electrical Engineering & Computer Sciences, University of California, Berkeley, CA, USA
   index: 1
 - name: Department of Genome Sciences, University of Washington, Seattle, WA, USA
   index: 2
 - name: Eccles Institute of Human Genetics, University of Utah, Salt Lake City, UT, USA
   index: 3
 - name: The Roslin Institute, Royal (Dick) School of Veterinary Studies, University of Edinburgh, Easter Bush Campus, Midlothian, UK
   index: 4
date: 18 January 2023
bibliography: paper.bib
---

# Summary

Genomics studies of the germline mutation process seek to elucidate the mutational forces that drive genetic variation and provide the raw material for adaptive evolution.
Germline mutations arise from spontaneous DNA damage caused by environmental mutagens, or errors in DNA replication.
Populations and species may experience distinct histories of mutational input due to differences in environmental exposure, life history (i.e. generation time), and heritable variation in the machinery controlling DNA replication fidelity.

Diverse mutational processes can reveal themselves in tell-tale mutation signatures left in the nucleotide sequence contexts in which they tend to be active.
Population genomics has given increasing attention to nucleotide sequence context in the study of the germline mutation process (reviewed in @Carlson2020-hb).
Single-nucleotide polymorphisms (SNPs) can be assigned to *mutation types* by the ancestral and derived nucleotide states and a window of local nucleotide context in the ancestral background.
The *mutation spectrum* of an individual or population is the relative distribution of these mutation types.

Inter- and intra-specific germline mutation spectrum variation has revealed a dynamic and evolving germline mutation process shaping modern genomic diversity.
Parsing mutation spectra temporally (i.e. according to allele frequency) and spatially (i.e. in different genomic compartments) has revealed the history and present of various mutational processes, and applying such analysis to *de novo* mutation data may be clinically informative for rare or undiagnosed genetic diseases.

Here we describe `mutyper`, a command-line utility and Python package that uses an ancestral genome estimate to assign mutation types to SNP data, compute mutation spectra for individuals and populations, and compute sample frequency spectra resolved by mutation type for population genetic inference.
Documentation is provided at [https://harrispopgen.github.io/mutyper](); source code is available at [https://github.com/harrispopgen/mutyper]().

# Statement of need

Despite many exciting findings to date, there is a lack of software for general-purpose germline mutation type partitioning and mutation spectrum generation from population-scale genomic variation data.
There is a need for efficient and well-tested software for both larger bioinformatics pipelines and simple exploratory analysis.
To address this need, we developed `mutyper`, an open-source command-line utility and Python package.

The literature on cancer somatic mutation signatures includes several software tools (many implemented as R packages) for clustering and dimensionality reduction that are not directly amenable to population-scale germline variation data, but the package `helmsman` [@Carlson2018-uq] enables interoperability with these tools.
Complementing this integrative work, `mutyper` is a flexible and extensible software package designed for population genomics researchers to generate the raw material needed to advance new analyses of germline mutation spectrum variation.

The literature has seen a recent and rapid expansion of studies of germline mutation spectrum variation and evolution.
This motivates the development of computational tools to make these analyses more standardized, reproducible, and accessible; `mutyper` meets this need.

# Implementation

## CLI

`mutyper` is a Python package with a command-line interface (CLI) that implements the core functionality of assigning ancestral mutation types to SNPs that are input (or piped) in VCF/BCF format.
Fast processing of VCF input [@Danecek2011-ng] is achieved with `cyvcf2` [@Pedersen2017-xu], and mutation types are assigned via the INFO field for each variant, e.g. by including a key-value pair such as `mutation_type=GAG>GTG`.
The user may specify the $k$-mer context size desired (e.g. 3 for triplet mutation types and spectra).
As in previous work, mutation types are, by default, collapsed by reverse complementation such that the ancestral state is either `A` or `C`.
Alternatively, a BED file can be supplied to define the strand orientation for nucleotide context at each site (e.g. based on direction of replication or transcription).

To polarize ancestral and derived allelic states, and define ancestral $k$-mer backgrounds, an input FASTA defining the ancestral genome estimate is required.
Mutyper uses the package `pyfaidx` [@Shirley2015-nf] for fast random access to ancestral genomic content, with minimal memory requirements.
Ancestral genomes can be specified by various means.
The ancestral FASTA sequence provided by the 1000 Genomes Project [@1000_Genomes_Project_Consortium2015-ek] was estimated from a multi-species alignment using `ortheus` [@Paten2008-ny].
In such a case, the ancestral FASTA can be passed to mutyper directly.
Alternatively, ancestral states can be simply estimated by polarizing SNPs using an outgroup genome aligned to the reference (e.g. the chimp genome liftover to the human genome).

More complex analyses can be achieved by filtering input SNPs (i.e. with `bcftools` [@Li2011-ca]), and piping to `mutyper`.
The `mutyper` command-line utility is fully compatible with standard command-line pipelines for filtering SNPs or samples, masking regions, and merging/concatenating VCFs.
In addition to this core functionality, the CLI includes several other subcommands that facilitate research that aims to characterize modern mutation spectrum variation, and infer its evolutionary history.

## Python API

The `mutyper` Python API enables one to perform the functions above in an interactive notebook session, or to implement custom analyses of mutation type data by interfacing with the strong ecosystem of scientific computing packages available in Python.
For example, dimensionality reduction (such as principal components analysis or non-negative matrix factorization) is often used to summarize mutation spectra, and the `scikit-learn` package [@scikit-learn] can be used in conjunction with the `mutyper` API for this purpose.
The `mutyper` API produces mutation spectra or SFS matrices as `pandas` data frames [@mckinney-proc-scipy-2010], which can be easily manipulated, visualized, and analyzed with standard python scientific computing packages.

# Applications

`mutyper` was first used by @DeWitt2021 alongside the Python package `mushi` to infer mutation rate histories from mutation spectra using coalescent theory.
@sasani2022 used `mutyper` in work reporting the discovery of a mutator allele in a unique mouse model system.
@Vollger2022 used `mutyper` to analyze long-read sequencing data from humans, finding elevated mutation rates and distinct mutation spectra in segmentally duplicated regions.
As of this writing, `mutyper` is being used in several ongoing studies in multiple labs.

# Acknowledgements

Jedidiah Carlson and Sarah Hilton provided helpful comments.
WSD was supported by the National Institute Of Allergy And Infectious Diseases (F31AI150163), and a Fellowship in Understanding Dynamic and Multi-scale Systems from the James S. McDonnell Foundation.
AT has been supported by the Institute Strategic Programme Grant BBS/E/D/10002070 from the Biotechnology and Biological Sciences Research Council (BBSRC).

# References