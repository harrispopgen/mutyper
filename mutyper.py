#! /usr/bin/env python

import cyvcf2
import argparse
import sys
from collections import defaultdict, Counter
import pandas as pd
import signal
import numpy as np

import ancestor


def setup_ancestor(args):
    """utility for initializing an Ancestor object for use in different '
    'subroutines"""
    return ancestor.Ancestor(args.fasta, k=args.k, target=args.target,
                             strand_file=args.strand_file,
                             key_function=lambda x:
                             x.split(args.sep)[args.chrom_pos],
                             read_ahead=10000,
                             sequence_always_upper=(not args.strict))


def variants(args):
    """subroutine for variants subcommand
    """
    ancestor = setup_ancestor(args)

    vcf = cyvcf2.VCF(args.vcf)
    vcf.add_info_to_header({'ID': 'mutation_type',
                            'Description': f'ancestral {args.k}-mer mutation '
                                           'type',
                            'Type': 'Character', 'Number': '1'})
    vcf_writer = cyvcf2.Writer('-', vcf)
    vcf_writer.write_header()
    for variant in vcf:
        # biallelic snps only
        if not (variant.is_snp and len(variant.ALT) == 1):
            continue
        # mutation type
        anc_kmer, der_kmer = ancestor.mutation_type(variant.CHROM,
                                                    variant.start, variant.REF,
                                                    variant.ALT[0])
        if anc_kmer is None or der_kmer is None:
            continue
        mutation_type = f'{anc_kmer}>{der_kmer}'
        variant.INFO['mutation_type'] = mutation_type
        # ancestral allele
        AA = ancestor[variant.CHROM][variant.start]
        # flip the alternative allele count if alternative allele is ancestral
        if variant.ALT[0] == AA:
            variant.INFO['AC'] = variant.INFO['AN'] - variant.INFO['AC']
        elif not variant.REF == AA:
            raise ValueError(f'ancestral allele {AA} is not equal to '
                             f'reference {variant.REF} or alternative '
                             f'{variant.ALT[0]}')
        # set REF to ancestral allele and ALT to derived allele
        variant.REF = anc_kmer[ancestor.target]
        variant.ALT = der_kmer[ancestor.target]
        vcf_writer.write_record(variant)
        # this line required to exit on a SIGTERM in a pipe, e.g. from head
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def targets(args):
    """subroutine for targets subcommand
    """
    ancestor = setup_ancestor(args)

    if args.bed == '-':
        args.bed = sys.stdin

    sizes = ancestor.targets(args.bed)

    try:
        for kmer in sorted(sizes):
            print(f'{kmer}\t{sizes[kmer]}')
    except BrokenPipeError:
        pass


def spectra(args):
    """subroutine for spectra subcommand
    """
    vcf = cyvcf2.VCF(args.vcf, gts012=True)

    spectra_data = defaultdict(lambda: np.zeros_like(vcf.samples, dtype=int))

    for variant in vcf:
        spectra_data[variant.INFO['mutation_type']] += variant.gt_types

    spectra = pd.DataFrame(spectra_data,
                           vcf.samples).reindex(sorted(spectra_data),
                                                axis='columns')
    try:
        print(spectra.to_csv(sep='\t', index=True,
                             index_label='sample'))
    except BrokenPipeError:
        pass


def ksfs(args):
    """subroutine for ksfs subcommand
    """
    vcf = cyvcf2.VCF(args.vcf)

    ksfs_data = defaultdict(lambda: Counter())
    AN = None
    for variant in vcf:
        # AN must be the same for all sites (no missing genotypes)
        if AN is not None and variant.INFO['AN'] != AN:
            raise ValueError(f'different AN {variant.INFO["AN"]} and {AN}'
                             ' indicates missing genotypes')
        AN = variant.INFO['AN']
        ksfs_data[variant.INFO['mutation_type']][variant.INFO['AC']] += 1

    index = range(1, AN)
    for mutation_type in sorted(ksfs_data):
        ksfs_data[mutation_type] = [ksfs_data[mutation_type][ac]
                                    for ac in index]
    ksfs = pd.DataFrame(ksfs_data, index).reindex(sorted(ksfs_data),
                                                  axis='columns')
    try:
        print(ksfs.to_csv(sep='\t', index=True,
                          index_label='sample_frequency'))
    except BrokenPipeError:
        pass


def main():
    """
    usage: python mutyper.py -h"""

    parser = argparse.ArgumentParser(
        description='mutyper: ancestral 𝑘-mer mutation types for variant data')
    subparsers = parser.add_subparsers(
        title='subcommands', description='specify one of these', required=True,
        help='additional help available for each subcommand')
    parser_variants = subparsers.add_parser(
        'variants', description='adds mutation_type to VCF/BCF INFO, polarizes'
                                ' REF/ALT/AC according to ancestral/derived '
                                'states, and stream to stdout')
    parser_targets = subparsers.add_parser(
        'targets', description='compute 𝑘-mer target sizes and stream to '
                               'stdout')
    parser_spectra = subparsers.add_parser(
        'spectra', description='compute mutation spectra for each sample in '
                               'VCF/BCF with mutation_type data and stream to'
                               ' stdout')
    parser_ksfs = subparsers.add_parser(
        'ksfs', description='compute sample frequency spectrum for each '
                            'mutation type from a VCF/BCF file with '
                            'mutation_type data (i.e. output from variants '
                            'subcommand ) and stream to stdout')

    # arguments that require FASTA input
    for sub_parser in (parser_variants, parser_targets):
        sub_parser.add_argument('fasta', type=str,
                                help='path to ancestral FASTA')
        sub_parser.add_argument('--k', type=int, default=3,
                                help='k-mer context size (default 3)')
        sub_parser.add_argument('--target', type=int, default=None,
                                help='0-based mutation target position in kmer'
                                     ' (default middle)')
        sub_parser.add_argument('--sep', type=str, default=':',
                                help='field delimiter in FASTA headers '
                                '(default ":")')
        sub_parser.add_argument('--chrom_pos', type=int, default=2,
                                help='0-based chromosome field position in '
                                'FASTA headers (default 2)')
        sub_parser.add_argument('--strand_file', type=str, default=None,
                                help='path to bed file with regions where '
                                     'reverse strand defines mutation '
                                     'context, e.g. direction of replication '
                                     'or transcription (default collapse '
                                     'reverse complements)')
        sub_parser.add_argument('--strict', action='store_true',
                                help='only uppercase nucleotides in FASTA '
                                     'considered ancestrally identified')

    # subcommands that require VCF input
    for sub_parser in (parser_variants, parser_spectra, parser_ksfs):
        sub_parser.add_argument('vcf', type=str,
                                help='VCF/BCF file created by the variants '
                                     'subcommand ("-" for stdin)')

    # arguments specific to variants subcommand
    parser_variants.set_defaults(func=variants)

    # arguments specific to targets subcommand
    parser_targets.add_argument('--bed', type=str, default=None,
                                help='path to BED file mask ("-" for stdin)')
    parser_targets.set_defaults(func=targets)

    # arguments specific to spectra subcommand
    parser_spectra.set_defaults(func=spectra)

    # arguments specific to ksfs subcommand
    parser_ksfs.set_defaults(func=ksfs)

    args = parser.parse_args()

    # execute the default subroutine for the chosen subcommand
    args.func(args)


if __name__ == '__main__':
    main()
