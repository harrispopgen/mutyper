mutyper
====

A Python3 package for annotating the local ancestral sequence context of biallelic SNPs, featuring

- A Python API
- A CLI for integration in pipelines analyzing VCF/BCF data


Dependencies
---

Dependences are listed in [env.yml](). You can set up a [conda](https://docs.conda.io/en/latest/) environment with
```bash
$ conda env create -f env.yml
```
and activate your new environment with
```bash
$ conda activate mutyper
```

Command line usage
---
### main help and subcommand list
```bash
$ python mutyper.py -h
usage: mutyper.py [-h] {variants,targets,spectra,ksfs} ...

mutyper: ancestral kmer mutation types for variant data

optional arguments:
  -h, --help            show this help message and exit

subcommands:
  specify one of these

  {variants,targets,spectra,ksfs}
                        additional help available for each subcommand
```

### `variants` subcommand
```bash
$ python mutyper.py variants -h
usage: mutyper.py variants [-h] [--k K] [--target TARGET] [--sep SEP]
                           [--chrom_pos CHROM_POS] [--strand_file STRAND_FILE]
                           [--strict]
                           fasta vcf

adds mutation_type to VCF/BCF INFO, polarizes REF/ALT/AC according to
ancestral/derived states, and stream to stdout

positional arguments:
  fasta                 path to ancestral FASTA
  vcf                   VCF/BCF file created by the variants subcommand ("-"
                        for stdin)

optional arguments:
  -h, --help            show this help message and exit
  --k K                 k-mer context size (default 3)
  --target TARGET       0-based mutation target position in kmer (default
                        middle)
  --sep SEP             field delimiter in FASTA headers (default ":")
  --chrom_pos CHROM_POS
                        0-based chromosome field position in FASTA headers
                        (default 2)
  --strand_file STRAND_FILE
                        path to bed file with regions where reverse strand
                        defines mutation context, e.g. direction of
                        replication or transcription (default collapse reverse
                        complements)
  --strict              only uppercase nucleotides in FASTA considered
                        ancestrally identified
```

### `targets` subcommand
```bash
$ python mutyper.py targets -h
usage: mutyper.py targets [-h] [--k K] [--target TARGET] [--sep SEP]
                          [--chrom_pos CHROM_POS] [--strand_file STRAND_FILE]
                          [--strict] [--bed BED]
                          fasta

compute 𝑘-mer target sizes and stream to stdout

positional arguments:
  fasta                 path to ancestral FASTA

optional arguments:
  -h, --help            show this help message and exit
  --k K                 k-mer context size (default 3)
  --target TARGET       0-based mutation target position in kmer (default
                        middle)
  --sep SEP             field delimiter in FASTA headers (default ":")
  --chrom_pos CHROM_POS
                        0-based chromosome field position in FASTA headers
                        (default 2)
  --strand_file STRAND_FILE
                        path to bed file with regions where reverse strand
                        defines mutation context, e.g. direction of
                        replication or transcription (default collapse reverse
                        complements)
  --strict              only uppercase nucleotides in FASTA considered
                        ancestrally identified
  --bed BED             path to BED file mask ("-" for stdin)

```

### `spectra` subcommand
```bash
$ python mutyper.py spectra -h
usage: mutyper.py spectra [-h] vcf

compute mutation spectra for each sample in VCF/BCF with mutation_type data
and stream to stdout

positional arguments:
  vcf         VCF/BCF file created by the variants subcommand ("-" for stdin)

optional arguments:
  -h, --help  show this help message and exit
```

### `ksfs` subcommand
```bash
$ python mutyper.py ksfs -h   
usage: mutyper.py ksfs [-h] vcf

compute sample frequency spectrum for each mutation type from a VCF/BCF file
with mutation_type data (i.e. output from variants subcommand ) and stream to
stdout

positional arguments:
  vcf         VCF/BCF file created by the variants subcommand ("-" for stdin)

optional arguments:
  -h, --help  show this help message and exit
```
