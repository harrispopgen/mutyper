#! /usr/bin/env python

import pyfaidx
import re
from Bio.Seq import reverse_complement
from typing import Generator, Tuple, Dict, Union, TextIO
from collections import Counter


class Ancestor(pyfaidx.Fasta):
    def __init__(self, fasta: str, k: int = 3, target: int = None,
                 strand_file: str = None, **kwargs):
        """ancestral state of a chromosome

        fasta: path to ancestral sequence FASTA
        k: the size of the context window (default 3)
        target: which position for the site within the kmer (default middle)
        strand_file: path to bed file with regions where reverse strand defines
                     mutation context, e.g. direction of replication or
                     transcription (default collapse reverse complements)
        kwargs: additional keyword arguments passed to base class. Useful
                ones are key_function (for chromosome name parsing), read_ahead
                (for buffering), and sequence_always_upper (to allow lowercase
                nucleotides to be considered ancestrally identified)
        """
        super(Ancestor, self).__init__(fasta, **kwargs)
        if target is None:
            assert k % 2 == 1, f'k = {k} must be odd for default middle target'
            target = k // 2
        else:
            raise NotImplementedError('target must be None (default)')
        assert 0 <= target < k
        self.target = target
        self.k = k
        if strand_file is None:
            self.strandedness = None
            if self.target != self.k // 2:
                raise ValueError(f'non-central target {self.target} requires '
                                 'strand_file')
        else:
            raise NotImplementedError('strand_file argument must be None')

    def mutation_type(self, chrom: str,
                      pos: int, ref: str, alt: str) -> Tuple[str, str]:
        """mutation type of a given snp, oriented or collapsed by strand
        returns a tuple of ancestral and derived kmers

        chrom: FASTA record chromosome identifier
        pos: position (0-based)
        ref: reference allele (A, C, G, or T)
        alt: alternative allele (A, C, G, or T)
        """
        # ancestral state
        anc = self[chrom][pos].seq
        # derived state
        if anc == ref:
            der = alt
        elif anc == alt:
            der = ref
        else:
            # infinite sites violation
            return None, None
        start = pos - self.target
        assert start >= 0
        end = pos + self.k - self.target
        assert start <= end

        context = self[chrom][start:end]
        anc_kmer = f'{context[:self.target]}{anc}{context[(self.target + 1):]}'
        der_kmer = f'{context[:self.target]}{der}{context[(self.target + 1):]}'

        if not re.match('^[ACGT]+$', anc_kmer) or not re.match('^[ACGT]+$',
                                                               der_kmer):
            return None, None

        if self.strandedness is None:
            if anc in 'AC':
                return anc_kmer, der_kmer
            elif anc in 'TG':
                return (reverse_complement(anc_kmer),
                        reverse_complement(der_kmer))
            else:
                raise ValueError('there is a bug if you got here')
        else:
            raise NotImplementedError('self.strandedness must be None')

    def region_contexts(self, chrom: str,
                        start: int = None,
                        end: int = None) -> Generator[str, None, None]:
        """ancestral context of each site in a BED style region (0-based,
        half-open), oriented according to self.strandedness or collapsed by
        reverse complementation (returns None if ancestral state at target not
        in capital ACGT)

        chrom: chromosome name
        start: region start position (default to chromsome start)
        end: region end position (default to chromsome end)
        """
        # NOTE: only valid for central target
        if start is None:
            start = self.target
        if end is None:
            end = len(self[chrom]) - self.target
        # we want to access the FASTA as few times as possible
        region_seq = self[chrom][(start - self.target):
                                 (end + self.k - self.target)]
        for i in range(end - start):
            context = region_seq[i:(i + self.k)].seq
            if not re.match('^[ACGT]+$', context):
                yield None
            elif self.strandedness is None:
                if context[self.target] in 'AC':
                    yield context
                elif context[self.target] in 'TG':
                    yield reverse_complement(context)
                else:
                    raise ValueError('there is a bug if you got here')
            else:
                raise NotImplementedError('self.strandedness must be None')

    def targets(self,
                bed: Union[str, TextIO] = None) -> Dict[str, int]:
        """return a dictionary of the number of sites of each k-mer

        bed: optional path to BED mask file, or I/O object"""
        sizes = Counter()
        if bed is None:
            for chrom in self.keys():
                sizes.update(self.region_contexts(chrom))
        else:
            if isinstance(bed, str):
                bed = open(bed, 'r')
            for line in bed:
                chrom, start, end = line.rstrip().split('\t')
                sizes.update(self.region_contexts(chrom, int(start), int(end)))
        del sizes[None]

        return sizes
