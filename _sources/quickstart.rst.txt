Quickstart
==========

This notebook demonstrates how to use ``mutyper`` for computing mutation
type data for VCF/BCF data

We will access some `example data <https://github.com/harrispopgen/mutyper/blob/master/docs/example_data/>`_ from a directory in the GitHub repo.

Python API demo
---------------

.. testsetup::

    # we move up to the root of the repo, so that example data paths in rendered docs
    # can be given wrt the root
    import os
    os.chdir('..')

.. testcode::

    import mutyper

Path to an ancestral FASTA file for human chromosome 1

.. testcode::

    fasta = 'docs/example_data/ancestor.fa'

Print the first 3 lines

.. testcode::

    with open(fasta) as f:
         for _ in range(3):
             print(f.readline().rstrip())

.. testoutput::

    >1
    TTTGTGGGAgACTATTCCTCCCATCTGCAACAGCTGCCCCTGCTGACTGCCCTTCTCTCCTCCCTCTCATCCCAGAGAAACAgGTCAGCTGGGAGCTTCT
    GCCCCCACTGCCTAGGGACCAACAGGGGCAGGAGGCAGTCACTGACCCCGAGACGTTTGCATCCTGCACAGCTAGAGATCCTTTATTAAAAGCACACTGT

The ``mutyper`` package has a single class ``Ancestor``. We instantiate
an ``Ancestor`` object using our FASTA file - We’re interested in 3-mer
context, which we can indicate with the ``k`` keyword argument (by
default it will be 3 though).

.. testcode::

    ancestor = mutyper.Ancestor(fasta, k=3)

We can inspect FASTA record names, showing that the FASTA contains two
records, chromosome ``'1'`` and ``'2'``.

.. testcode::

    print(ancestor.keys())

.. testoutput::

    odict_keys(['1', '2'])


and we can Pythonically slice sequences via fast random access without
loading into memory (leaning on
`pyfaidx <https://pythonhosted.org/pyfaidx/>`__ under the hood):

.. testcode::

    start = 100
    end = 300
    print(repr(ancestor['1'][start:end]))

.. testoutput::

    >1:101-300
    GCCCCCACTGCCTAGGGACCAACAGGGGCAGGAGGCAGTCACTGACCCCGAGACGTTTGCATCCTGCACAGCTAGAGATCCTTTATTAAAAGCACACTGTTGGTTTCTGCTCAGTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGgAG


We can also access these ``FastaRecord`` slices as string-like biopython
``Seq`` objects

.. testcode::

    print(ancestor['1'][start:end].seq)

.. testoutput::

    GCCCCCACTGCCTAGGGACCAACAGGGGCAGGAGGCAGTCACTGACCCCGAGACGTTTGCATCCTGCACAGCTAGAGATCCTTTATTAAAAGCACACTGTTGGTTTCTGCTCAGTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCATGgAG



The ``mutation_type`` method allows us to specify a SNP by the usual
CHROM, POS, REF, ALT, and returns the correctly polarized mutation type
as a tuple of ancestral kmer and derived kmer.

First, consider the triplet context at site 100 on chromosome 1:

.. testcode::

    print(repr(ancestor['1'][98:101]))

.. testoutput::

    >1:99-101
    CTG

Suppose a biallelic SNP at this site segregates in a population with
reference allele T and alternative allele A. The mutation type of this
SNP is

.. testcode::

    print(ancestor.mutation_type('1', 99, 'T', 'A'))

.. testoutput::

    ('CAG', 'CTG')


Note that, by default, mutation types are collapsed by reverse
complementation such that the ancestral state at the target site is A or
C. So in the above, the ancestral state T caused the ancestral triplet
CTG to be reversed and complemented to CAG.

If the reference and alternative states are both distinct from the
ancestral state, the site is not biallelic (an infinite sites
violation), and the mutation type is returned as ambiguous.

.. testcode::

    print(ancestor.mutation_type('1', 99, 'C', 'A'))

.. testoutput::

    (None, None)

The ``region_contexts`` method returns a generator of ancestral contexts
(by default collapsed as described above) over the positions in a region
specified as in BED file format CHROM, START, END

.. testcode::

    start = 50
    end = 100
    for context in ancestor.region_contexts('1', start, end):
        print(context, end=' ')

.. testoutput::
    :options: +NORMALIZE_WHITESPACE

    CCC CCT AAG GAA TCT GAG TCT GAG TCC CCT GAG TCC CCC CCT GAG TCT GAG TCA CAT GAT TCC CCC CCA CAG TCT GAG TCT GAA AAA AAC ACA None None None GAC TCA CAG GCT GCT CAG CCA CCC TCC GAG GCT GCT AAG GAA TCT CAG 

Note that FASTA sites that are not nucleotide characters have
non-identified ancestral state, so context is ``None``.

This generator method is used by the the ``targets`` method to compute
the masked genomic target size for each :math:`k`-mer context from a BED
mask file. This may be useful for normalizing spectra in different
genomic regions, or calibrating mutation rates.

.. testcode::

    bed = 'docs/example_data/mask.bed'

.. testcode::

    with open(bed) as f:
         for _ in range(3):
             print(f.readline().rstrip())

.. testoutput::
    :options: +NORMALIZE_WHITESPACE

    1	25	50
    1	150	175
    2	50	100


Under the hood it loops over BED file entries and updates a ``Counter``
object for the different triplets

.. testcode::

    print(ancestor.targets(bed))

.. testoutput::

    Counter({'CAG': 10, 'GCA': 9, 'GCT': 7, 'TCT': 7, 'AAA': 7, 'CAT': 5, 'GAT': 5, 'GCC': 4, 'GAG': 4, 'TAG': 4, 'TAA': 4, 'CCA': 4, 'AAC': 3, 'CCC': 3, 'CCT': 3, 'TCC': 3, 'AAT': 3, 'CAA': 2, 'ACA': 2, 'TCA': 2, 'GAC': 2, 'ACG': 2, 'ACT': 1, 'CAC': 1, 'GAA': 1, 'AAG': 1, 'ACC': 1})



We can also call it without the BED file mask to compute the unmasked
target sizes of the ancestral FASTA

.. testcode::

    print(ancestor.targets())

.. testoutput::

    Counter({'AAA': 6197, 'CAG': 4423, 'TCT': 4319, 'AAT': 4240, 'ACA': 3920, 'GAA': 3861, 'AAG': 3860, 'CCT': 3811, 'CCA': 3764, 'TCA': 3730, 'GAG': 3601, 'TAA': 3527, 'CAT': 3456, 'CAA': 3408, 'TCC': 3317, 'TAT': 3279, 'CAC': 3260, 'CCC': 3234, 'GCA': 3182, 'ACT': 2998, 'GCT': 2972, 'GCC': 2694, 'AAC': 2688, 'ACC': 2378, 'GAT': 2370, 'TAG': 2210, 'GAC': 1998, 'TAC': 1927, 'CCG': 717, 'GCG': 592, 'ACG': 554, 'TCG': 450})



Command line interface demo
---------------------------

Display usage information for the mutyper command

.. code:: console

    $ mutyper -h


.. parsed-literal::

    usage: mutyper [-h] {ancestor,variants,targets,spectra,ksfs} ...
    
    mutyper: ancestral kmer mutation types for variant data
    
    options:
      -h, --help            show this help message and exit
    
    subcommands:
      specify one of these
    
      {ancestor,variants,targets,spectra,ksfs}
                            additional help available for each subcommand


``variants`` subcommand
~~~~~~~~~~~~~~~~~~~~~~~

Usage:

.. code:: console

    $ mutyper variants -h


.. parsed-literal::

    usage: mutyper variants [-h] [--verbose] [--k K] [--target TARGET] [--sep SEP]
                            [--chrom_pos CHROM_POS] [--strand_file STRAND_FILE]
                            [--strict]
                            fasta vcf
    
    adds mutation_type to VCF/BCF INFO, polarizes REF/ALT/AC according to
    ancestral/derived states, and stream to stdout
    
    positional arguments:
      fasta                 path to ancestral FASTA
      vcf                   VCF/BCF file, usually for a single chromosome ("-" for
                            stdin)
    
    options:
      -h, --help            show this help message and exit
      --verbose             increase logging verbosity
      --k K                 k-mer context size (default 3)
      --target TARGET       0-based mutation target position in kmer (default
                            middle)
      --sep SEP             field delimiter in FASTA headers (default ":")
      --chrom_pos CHROM_POS
                            0-based chromosome field position in FASTA headers
                            (default 0)
      --strand_file STRAND_FILE
                            path to bed file with regions where reverse strand
                            defines mutation context, e.g. direction of
                            replication or transcription (default collapse reverse
                            complements)
      --strict              only uppercase nucleotides in FASTA considered
                            ancestrally identified


Path to FASTA and a truncated VCF file for chromosome 1 from the 1000 Genomes
Project.

.. code:: console

    $ fasta = 'docs/example_data/ancestral.fa'
    $ vcf = 'docs/example_data/snps.vcf'

We'll use
`bcftools <http://samtools.github.io/bcftools/bcftools.html#view>`__ to
display the first few variants, omitting header and genotype
information.

.. code:: console

    $ bcftools view -HG $vcf | head -5


.. parsed-literal::

    1       10177   .       A       AC      100     PASS    AC=2130;AF=0.425319;AN=5008;NS=2504;DP=103152;EAS_AF=0.3363;AMR_AF=0.3602;AFR_AF=0.4909;EUR_AF=0.4056;SAS_AF=0.4949;AA=|||unknown(NO_COVERAGE)
    1       10235   .       T       TA      100     PASS    AC=6;AF=0.00119808;AN=5008;NS=2504;DP=78015;EAS_AF=0;AMR_AF=0.0014;AFR_AF=0;EUR_AF=0;SAS_AF=0.0051;AA=|||unknown(NO_COVERAGE)
    1       10352   rs145072688     T       TA      100     PASS    AC=2191;AF=0.4375;AN=5008;NS=2504;DP=88915;EAS_AF=0.4306;AMR_AF=0.4107;AFR_AF=0.4788;EUR_AF=0.4264;SAS_AF=0.4192;AA=|||unknown(NO_COVERAGE)
    1       10505   .       A       T       100     PASS    AC=1;AF=0.000199681;AN=5008;NS=2504;DP=9632;EAS_AF=0;AMR_AF=0;AFR_AF=0.0008;EUR_AF=0;SAS_AF=0;AA=.|||
    1       10506   .       C       G       100     PASS    AC=1;AF=0.000199681;AN=5008;NS=2504;DP=9676;EAS_AF=0;AMR_AF=0;AFR_AF=0.0008;EUR_AF=0;SAS_AF=0;AA=.|||


Now we use the ``variants`` subcommand, and pipe to ``bcftools`` for
display. Notice there is an additional INFO field ``mutation_type``, and
only biallelic SNPs are output. Also note that REF/ALT state and INFO/AC
fields are polarized according to ancestral/derived allele.

.. code:: console

    $ mutyper variants $fasta $vcf | bcftools view -HG | head -5


.. parsed-literal::

    1       10505   .       A       T       100     PASS    AC=5007;AF=0.9998;AN=5008;NS=2504;DP=9632;EAS_AF=0;AMR_AF=0;AFR_AF=0.0008;EUR_AF=0;SAS_AF=0;AA=.|||;mutation_type=GAG>GTG
    1       10506   .       C       G       100     PASS    AC=1;AF=0.000199681;AN=5008;NS=2504;DP=9676;EAS_AF=0;AMR_AF=0;AFR_AF=0.0008;EUR_AF=0;SAS_AF=0;AA=.|||;mutation_type=TCT>TGT
    1       10542   .       A       G       100     PASS    AC=5007;AF=0.9998;AN=5008;NS=2504;DP=9007;EAS_AF=0.001;AMR_AF=0;AFR_AF=0;EUR_AF=0;SAS_AF=0;AA=.|||;mutation_type=CAG>CGG
    1       10642   .       A       G       100     PASS    AC=4987;AF=0.995807;AN=5008;NS=2504;DP=1360;EAS_AF=0.003;AMR_AF=0.0014;AFR_AF=0.0129;EUR_AF=0;SAS_AF=0;AA=.|||;mutation_type=CAT>CGT
    1       11008   .       C       G       100     PASS    AC=441;AF=0.0880591;AN=5008;NS=2504;DP=2232;EAS_AF=0.0367;AMR_AF=0.0965;AFR_AF=0.1346;EUR_AF=0.0885;SAS_AF=0.0716;AA=.|||;mutation_type=CCT>CGT

``targets`` subcommand
~~~~~~~~~~~~~~~~~~~~~~

Usage:

.. code:: console

    $ mutyper targets -h


.. parsed-literal::

    usage: mutyper targets [-h] [--verbose] [--k K] [--target TARGET] [--sep SEP]
                           [--chrom_pos CHROM_POS] [--strand_file STRAND_FILE]
                           [--strict] [--bed BED]
                           fasta
    
    compute 𝑘-mer target sizes and stream to stdout
    
    positional arguments:
      fasta                 path to ancestral FASTA
    
    options:
      -h, --help            show this help message and exit
      --verbose             increase logging verbosity
      --k K                 k-mer context size (default 3)
      --target TARGET       0-based mutation target position in kmer (default
                            middle)
      --sep SEP             field delimiter in FASTA headers (default ":")
      --chrom_pos CHROM_POS
                            0-based chromosome field position in FASTA headers
                            (default 0)
      --strand_file STRAND_FILE
                            path to bed file with regions where reverse strand
                            defines mutation context, e.g. direction of
                            replication or transcription (default collapse reverse
                            complements)
      --strict              only uppercase nucleotides in FASTA considered
                            ancestrally identified
      --bed BED             path to BED file mask ("-" for stdin)


.. code:: console

    $ mutyper targets $fasta | head


.. parsed-literal::

    AAA	6331
    AAC	2775
    AAG	3947
    AAT	4335
    ACA	4066
    ACC	2440
    ACG	632
    ACT	3082
    CAA	3506
    CAC	3391


The optional argument ``--bed`` resticts target size computation based
on a BED mask file, and is useful of normalizing mutation spectra
according to genomic region, or calibrating mutaiton rates.

.. code:: console

    $ mutyper targets $fasta --bed $bed | head


.. parsed-literal::

    AAA	7
    AAC	3
    AAG	1
    AAT	3
    ACA	2
    ACC	1
    ACG	2
    ACT	1
    CAA	2
    CAC	1


``spectra`` subcommand
~~~~~~~~~~~~~~~~~~~~~~

Usage:

.. code:: console

    $ mutyper spectra -h


.. parsed-literal::

    usage: mutyper spectra [-h] [--verbose] [--population] [--randomize] vcf
    
    compute mutation spectra for each sample in VCF/BCF with mutation_type data
    and stream to stdout
    
    positional arguments:
      vcf           VCF/BCF file, usually for a single chromosome ("-" for stdin)
    
    options:
      -h, --help    show this help message and exit
      --verbose     increase logging verbosity
      --population  population-wise spectrum, instead of individual
      --randomize   randomly assign mutation to a single haplotype


We first add mutation type information with the ``variants`` subcommand,
then pipe into the ``spectra`` subcommand to generate mutation spectra
for all samples, and finally crop the output with ``head`` and ``cut``
for display.

.. code:: console

    $ mutyper variants $fasta $vcf | mutyper spectra - | head | cut -f-5


.. parsed-literal::

    sample	AAA>AGA	AAT>AGT	ACA>AAA	ACA>AGA
    HG00096	2	2	2	2
    HG00097	2	2	2	2
    HG00099	2	2	2	2
    HG00100	2	2	2	2
    HG00101	2	2	2	2
    HG00102	2	2	2	2
    HG00103	2	2	2	2
    HG00105	2	2	2	2
    HG00106	2	2	2	2


Use the ``--population`` to get the spectrum for whole population,
rather than each individual

.. code:: console

    $ mutyper variants $fasta $vcf | mutyper spectra - --population | head | cut -f-5


.. parsed-literal::

    AAA>AGA	AAT>AGT	ACA>AAA	ACA>AGA	ACA>ATA
    2	1	1	2	3
    


``ksfs`` subcommand
~~~~~~~~~~~~~~~~~~~

Usage:

.. code:: console

    $ mutyper ksfs -h


.. parsed-literal::

    usage: mutyper ksfs [-h] [--verbose] vcf
    
    compute sample frequency spectrum for each mutation type from a VCF/BCF file
    with mutation_type data (i.e. output from variants subcommand ) and stream to
    stdout
    
    positional arguments:
      vcf         VCF/BCF file, usually for a single chromosome ("-" for stdin)
    
    options:
      -h, --help  show this help message and exit
      --verbose   increase logging verbosity


Similar to the previous subcommand, but now we are generating a sample
frequency spectrum (SFS) for each mutation type.

.. code:: console

    $ mutyper variants $fasta $vcf | mutyper ksfs - | head | cut -f-5


.. parsed-literal::

    sample_frequency	AAA>AGA	AAT>AGT	ACA>AAA	ACA>AGA
    1	0	0	0	0
    2	0	0	0	0
    3	0	0	0	0
    4	0	0	0	0
    5	0	0	0	0
    6	0	0	0	0
    7	0	0	0	0
    8	0	0	0	0
    9	0	0	0	0

