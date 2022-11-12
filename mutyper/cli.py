#! /usr/bin/env python

import cyvcf2
import argparse
import sys
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from shutil import copyfile, copyfileobj
import pyfaidx
from pyliftover import LiftOver
from Bio.Seq import reverse_complement
import logging
import gzip
from mutyper import ancestor


def setup_ancestor(args):
    """Utility for initializing an Ancestor object for use in different
    subroutines."""
    return ancestor.Ancestor(
        args.fasta,
        k=args.k,
        target=args.target,
        strand_file=args.strand_file,
        key_function=lambda x: x.split(args.sep)[args.chrom_pos],
        read_ahead=10000,
        sequence_always_upper=(not args.strict),
    )


def is_compressed(file):
    """Returns ``True`` if file is compressed."""
    f = open(file, "rb")
    # The first two bytes of a gzip file are: 1f 8b
    compressed = f.read(2) == b"\x1f\x8b"
    f.close()
    logging.info(f"{file} is {'compressed' if compressed else 'not compressed'}")
    return compressed


def copy_fasta(file, outfile):
    """Copy FASTA file to ``outfile`` and decompress if needed.

    This is necessary because pyfaidx does not support mutable
    compressed, so we cannot just copy the input if compressed.
    """
    # outfile cannot be compressed
    if outfile.endswith(".gz"):
        raise ValueError(
            f"{outfile} target cannot be compressed because pafaidx needs a mutable file type (remove .gz)!"
        )
    # check if input is compressed and make decompressed copy
    if is_compressed(file):
        with gzip.open(file, "r") as f_in, open(outfile, "wb") as f_out:
            copyfileobj(f_in, f_out)
    else:
        copyfile(file, outfile)


def ancestral_fasta(args):
    """Subroutine for ancestor subcommand."""
    # single chromosome fasta file for reference genome
    ref = pyfaidx.Fasta(args.reference, read_ahead=10000)
    # make a copy to build our ancestor for this chromosome
    copy_fasta(args.reference, args.output)
    anc = pyfaidx.Fasta(args.output, read_ahead=10000, mutable=True)
    # reference genome for outgroup species (all chromosomes)
    out = pyfaidx.Fasta(args.outgroup, read_ahead=10000)
    # outgroup to reference alignment chain file
    lo = LiftOver(args.chain)
    # snps database for the same chromosome
    vcf = cyvcf2.VCF(args.vcf)

    # change regions outside of callability mask to all N bases
    if args.bed:
        if args.bed == "-":
            bed = sys.stdin
        else:
            bed = open(args.bed, "r")
        last_end = 0
        for line in bed:
            chrom, start, end = line.rstrip().split("\t")[:3]
            start = int(start)
            anc[chrom][last_end:start] = "N" * (start - last_end)
            last_end = int(end)
        anc[chrom][last_end : len(anc[chrom])] = "N" * (len(anc[chrom]) - last_end)

    # variables for counting liftover cases
    num_vars = 0
    num_changed = 0
    num_unchanged = 0
    num_unknown = 0

    for variant in vcf:
        num_vars += 1
        # change variants that are not biallelic SNPs to N bases
        if not (variant.is_snp and len(variant.ALT) == 1):
            anc[variant.CHROM][variant.start : variant.end] = "N" * (
                variant.end - variant.start
            )
            num_unknown += 1
        else:
            out_coords = lo.convert_coordinate(variant.CHROM, variant.start)
            # change ambiguously aligning sites to N bases
            if out_coords is None or len(out_coords) != 1:
                anc[variant.CHROM][variant.start] = "N"
                num_unknown += 1
            else:
                if variant.REF != ref[variant.CHROM][variant.start].seq.upper():
                    raise ValueError(
                        f"variant reference allele {variant.REF} "
                        f"mismatches reference sequence "
                        f"{ref[variant.CHROM][variant.start]}"
                    )
                out_chromosome, out_position, out_strand = out_coords[0][:3]
                out_allele = out[out_chromosome][out_position].seq
                # if negative strand, take reverse complement base
                if out_strand == "-":
                    out_allele = reverse_complement(out_allele)
                # and finally, polarize
                if out_allele.upper() == variant.ALT[0]:
                    anc[variant.CHROM][variant.start] = out_allele
                    num_changed += 1
                elif out_allele.upper() != variant.REF:
                    # triallelic
                    anc[variant.CHROM][variant.start] = "N"
                    num_unknown += 1
                else:
                    num_unchanged += 1

    logging.info(f"{num_vars} total variant positions.")
    logging.info(f"{num_unchanged} variant positions unchanged.")
    logging.info(f"{num_changed} variant positions changed.")
    logging.info(
        f"{num_unknown} variant positions with missing or ambiguous liftOvers."
    )
    if num_unknown / num_vars > 0.5:
        logging.warning(
            f"{num_unknown / num_vars * 100:.2f}%) ({num_unknown}/{num_vars}) "
            f"of variant positions are missing or ambiguous liftOvers."
            f"Check that your chain file is correct."
        )


def variants(args):
    """Subroutine for variants subcommand."""
    ancestor = setup_ancestor(args)

    vcf = cyvcf2.VCF(args.vcf)
    vcf.add_info_to_header(
        {
            "ID": "mutation_type",
            "Description": f"ancestral {args.k}-mer mutation type",
            "Type": "Character",
            "Number": "A",
        }
    )
    print(vcf.raw_header, end="")
    num_vars = 0
    valid_ploidy = set([1, 2])
    for variant in vcf:
        num_vars += 1
        # biallelic snps only
        if not (variant.is_snp and len(variant.ALT) == 1):
            continue

        if variant.ploidy not in valid_ploidy:
            raise ValueError(
                f"invalid ploidy {variant.ploidy}, diploids and haploids only"
            )

        # mutation type as ancestral kmer and derived kmer
        anc_kmer, der_kmer = ancestor.mutation_type(
            variant.CHROM, variant.start, variant.REF, variant.ALT[0]
        )
        if anc_kmer is None or der_kmer is None:
            continue
        mutation_type = f"{anc_kmer}>{der_kmer}"
        variant.INFO["mutation_type"] = mutation_type

        # checks that all genotype elements are from the set of {-1,0,1}
        # each element in the last column is a 0,1 indicator for phasing status
        genotype_array = variant.genotype.array()
        unique_gts = set(np.unique(genotype_array[:, :-1]))
        if not unique_gts <= set([-1, 0, 1]):
            raise ValueError(f"invalid genotypes {unique_gts - set([-1,0,1])}")

        # ancestral allele
        AA = ancestor[variant.CHROM][variant.start].seq

        # polarize genotypes (and associated INFO) if alternative allele is
        # ancestral
        if variant.ALT[0] == AA:
            variant.INFO["AC"] = variant.INFO["AN"] - variant.INFO["AC"]
            variant.INFO["AF"] = variant.INFO["AC"] / variant.INFO["AN"]

            # cyvcf2 docs say we need to reassign genotypes like this for the
            # change to propagate (can't just update indexwise)
            genotype_array[:, :-1] = np.select(
                [genotype_array[:, :-1] == state for state in [-1, 0, 1]], [-1, 1, 0]
            )
            variant.genotypes = genotype_array
        elif not variant.REF == AA:
            raise ValueError(
                f"ancestral allele {AA} is not equal to "
                f"reference {variant.REF} or alternative "
                f"{variant.ALT[0]}"
            )
        # set REF to ancestral allele and ALT to derived allele
        variant.REF = anc_kmer[ancestor.target]
        variant.ALT = der_kmer[ancestor.target]
        print(str(variant), end="")

    if num_vars == 0:
        logging.warning("No variants processed. Check that input vcf is not empty.")


def targets(args):
    """Subroutine for targets subcommand."""
    ancestor = setup_ancestor(args)

    if args.bed == "-":
        args.bed = sys.stdin

    sizes = ancestor.targets(args.bed)

    try:
        for kmer in sorted(sizes):
            print(f"{kmer}\t{sizes[kmer]}")
    except BrokenPipeError:
        pass


def spectra(args):
    """Subroutine for spectra subcommand."""

    vcf = cyvcf2.VCF(args.vcf, strict_gt=True)

    def iterate_with_ambiguity_warning():
        """In several places we want to check for genotype ambiguity as we
        iterate over vcf variants, so we define this generator wrapper."""
        for variant in vcf:
            yield variant
            if variant.num_unknown:
                logging.warning(
                    "Ambiguous genotypes found! Continuing by assuming reference genotypes for these variants."
                )
                break
        yield from vcf

    if args.population:
        spectra_data = Counter()

        for variant in iterate_with_ambiguity_warning():
            if variant.aaf:
                spectra_data[variant.INFO["mutation_type"]] += 1

        spectra = pd.DataFrame(spectra_data, ["population"]).reindex(
            sorted(spectra_data), axis="columns"
        )
        try:
            print(spectra.to_csv(sep="\t", index=False))
        except BrokenPipeError:
            pass

    else:
        spectra_data = defaultdict(lambda: np.zeros_like(vcf.samples, dtype=int))
        if args.randomize:
            for variant in iterate_with_ambiguity_warning():
                counts = (variant.genotype.array()[:, :-1] == 1).sum(axis=1)
                rng = np.random.default_rng()
                random_haplotype = rng.choice(len(counts), p=counts / counts.sum())
                spectra_data[variant.INFO["mutation_type"]][random_haplotype] += 1.0
        else:
            for variant in iterate_with_ambiguity_warning():
                counts = (variant.genotype.array()[:, :-1] == 1).sum(axis=1)
                spectra_data[variant.INFO["mutation_type"]] += counts

        spectra = pd.DataFrame(spectra_data, vcf.samples).reindex(
            sorted(spectra_data), axis="columns"
        )
        try:
            print(spectra.to_csv(sep="\t", index=True, index_label="sample"))
        except BrokenPipeError:
            pass


def ksfs(args):
    """Subroutine for ksfs subcommand."""
    vcf = cyvcf2.VCF(args.vcf)

    ksfs_data = defaultdict(lambda: Counter())
    AN = None
    for variant in vcf:
        # AN must be the same for all sites (no missing genotypes)
        if AN is not None and variant.INFO["AN"] != AN:
            raise ValueError(
                f'different AN {variant.INFO["AN"]} and {AN}'
                " indicates missing genotypes"
            )
        AN = variant.INFO["AN"]
        ksfs_data[variant.INFO["mutation_type"]][variant.INFO["AC"]] += 1

    # exclude fixed sites AC=0, AC=AN
    index = range(1, AN)
    for mutation_type in sorted(ksfs_data):
        ksfs_data[mutation_type] = [ksfs_data[mutation_type][ac] for ac in index]
    ksfs = pd.DataFrame(ksfs_data, index).reindex(sorted(ksfs_data), axis="columns")
    try:
        print(ksfs.to_csv(sep="\t", index=True, index_label="sample_frequency"))
    except BrokenPipeError:
        pass


def get_parser():
    parser = argparse.ArgumentParser(
        description="mutyper: ancestral kmer mutation types for variant data"
    )
    subparsers = parser.add_subparsers(
        title="subcommands",
        description="specify one of these",
        required=True,
        help="additional help available for each subcommand",
    )
    parser_ancestor = subparsers.add_parser(
        "ancestor",
        description="create an ancestral FASTA file, using an "
        "outgroup alignment, and a database of SNPs "
        "with a callability mask",
    )
    parser_variants = subparsers.add_parser(
        "variants",
        description="adds mutation_type to VCF/BCF INFO, polarizes"
        " REF/ALT/AC according to ancestral/derived "
        "states, and stream to stdout",
    )
    parser_targets = subparsers.add_parser(
        "targets", description="compute 𝑘-mer target sizes and stream to stdout"
    )
    parser_spectra = subparsers.add_parser(
        "spectra",
        description="compute mutation spectra for each sample in "
        "VCF/BCF with mutation_type data and stream to"
        " stdout",
    )
    parser_ksfs = subparsers.add_parser(
        "ksfs",
        description="compute sample frequency spectrum for each "
        "mutation type from a VCF/BCF file with "
        "mutation_type data (i.e. output from variants "
        "subcommand ) and stream to stdout",
    )

    # arguments for all subparsers
    for sub_parser in (
        parser_ancestor,
        parser_variants,
        parser_targets,
        parser_spectra,
        parser_ksfs,
    ):
        sub_parser.add_argument(
            "--verbose", help="increase logging verbosity", action="store_true"
        )

    # arguments that require FASTA input
    for sub_parser in (parser_variants, parser_targets):
        sub_parser.add_argument("fasta", type=str, help="path to ancestral FASTA")
        sub_parser.add_argument(
            "--k", type=int, default=3, help="k-mer context size (default 3)"
        )
        sub_parser.add_argument(
            "--target",
            type=int,
            default=None,
            help="0-based mutation target position in kmer (default middle)",
        )
        sub_parser.add_argument(
            "--sep",
            type=str,
            default=":",
            help='field delimiter in FASTA headers (default ":")',
        )
        sub_parser.add_argument(
            "--chrom_pos",
            type=int,
            default=0,
            help="0-based chromosome field position in FASTA headers (default 0)",
        )
        sub_parser.add_argument(
            "--strand_file",
            type=str,
            default=None,
            help="path to bed file with regions where "
            "reverse strand defines mutation "
            "context, e.g. direction of replication "
            "or transcription (default collapse "
            "reverse complements)",
        )
        sub_parser.add_argument(
            "--strict",
            action="store_true",
            help="only uppercase nucleotides in FASTA "
            "considered ancestrally identified",
        )

    # subcommands that require VCF input
    for sub_parser in (parser_ancestor, parser_variants, parser_spectra, parser_ksfs):
        sub_parser.add_argument(
            "vcf",
            type=str,
            help="VCF/BCF file, usually for a single " 'chromosome ("-" for stdin)',
        )

    # subcommands that take BED input
    for sub_parser in (parser_ancestor, parser_targets):
        sub_parser.add_argument(
            "--bed",
            type=str,
            default=None,
            help='path to BED file mask ("-" for stdin)',
        )

    # arguments specific to ancestor subcommand
    parser_ancestor.add_argument(
        "reference", type=str, help="path to reference FASTA for one chromosome"
    )
    parser_ancestor.add_argument(
        "outgroup", type=str, help="path to outgroup genome FASTA"
    )
    parser_ancestor.add_argument(
        "chain",
        type=str,
        help="path to alignment chain file (reference to outgroup)",
    )
    parser_ancestor.add_argument(
        "output",
        type=str,
        help="path for output ancestral FASTA for this chromosome",
    )
    parser_ancestor.set_defaults(func=ancestral_fasta)

    # arguments specific to variants subcommand
    parser_variants.set_defaults(func=variants)

    # arguments specific to targets subcommand
    parser_targets.set_defaults(func=targets)

    # arguments specific to spectra subcommand
    parser_spectra.add_argument(
        "--population",
        action="store_true",
        help="population-wise spectrum, instead of individual",
    )
    parser_spectra.add_argument(
        "--randomize",
        action="store_true",
        help="randomly assign mutation to a single haplotype",
    )
    parser_spectra.set_defaults(func=spectra)

    # arguments specific to ksfs subcommand
    parser_ksfs.set_defaults(func=ksfs)

    return parser


def main(arg_list=None):
    args = get_parser().parse_args(arg_list)
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(format=log_format, level=log_level)
    args.func(args)
