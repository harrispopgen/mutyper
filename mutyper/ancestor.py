#! /usr/bin/env python

import pyfaidx
import re
from Bio.Seq import reverse_complement
from typing import Generator, Tuple, Dict, Union, TextIO
from collections import Counter, defaultdict
import bisect


class Ancestor(pyfaidx.Fasta):
    r"""ancestral state of a chromosome

    Args:
        fasta: path to ancestral sequence FASTA
        k: the size of the context window (default 3)
        target: which position for the site within the kmer (default middle)
        strand_file: path to bed file (or I/O object) with regions where
                     reverse strand defines mutation context, e.g. direction of
                     replication or transcription. Sites not in theses regions
                     are assigned forward strand context. If not provided,
                     collapse by reverse complement.
        kwargs: additional keyword arguments passed to base class. Useful
                ones are ``key_function`` (for chromosome name parsing),
                ``read_ahead`` (for buffering), and ``sequence_always_upper``
                (to allow lowercase nucleotides to be considered ancestrally
                identified)
        """
    acgt = set('ACGT')
    ac = set('AC')

    def __init__(self, fasta: str, k: int = 3, target: int = None,
                 strand_file: Union[str, TextIO] = None, **kwargs):
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
            self._revcomp_func = self._AC
            if self.target != self.k // 2:
                raise ValueError(f'non-central target {self.target} requires '
                                 'strand_file')
        else:
            self.strandedness = defaultdict(list)
            if isinstance(strand_file, str):
                bed = open(strand_file, 'r')
            for line in bed:
                chrom, start, end = line.rstrip().split('\t')
                bisect.insort(self.strandedness[chrom], (int(start), int(end)))
            self._revcomp_func = self._reverse_strand

    def _reverse_strand(self, chrom: str, pos: int):
        r"""return True if strand_file indicates reverse complementation at
        this site"""
        closest_idx = bisect.bisect(self.strandedness[chrom], (pos, -1))
        if pos < self.strandedness[chrom][closest_idx][1]:
            return True
        return False

    def _AC(self, chrom: str, pos: int):
        r"""return True if reverse complementation is needed at this site
        to get state A or C"""
        if self[chrom][pos].seq not in self.ac:
            return True
        return False

    def mutation_type(self, chrom: str,
                      pos: int, ref: str, alt: str) -> Tuple[str, str]:
        r"""mutation type of a given snp, oriented or collapsed by strand,
        returns a tuple of ancestral and derived kmers

        Args:
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

        if not self._revcomp_func(chrom, pos):
            return anc_kmer, der_kmer
        else:
            return (reverse_complement(anc_kmer),
                    reverse_complement(der_kmer))

    def region_contexts(self, chrom: str,
                        start: int = None,
                        end: int = None) -> Generator[str, None, None]:
        r"""ancestral context of each site in a BED style region (0-based,
        half-open), oriented according to self.strandedness or collapsed by
        reverse complementation (returns None if ancestral state at target not
        in capital ACGT)

        Args:
            chrom: chromosome name
            start: region start position (default to chromsome start)
            end: region end position (default to chromsome end)
        """
        if start is None:
            start = 0
        if end is None:
            end = len(self[chrom])
        for pos in range(start, end):
            if self[chrom][pos].seq not in self.acgt:
                yield None
                continue
            if not self._revcomp_func(chrom, pos):
                context_start = pos - self.target
                context_end = pos + self.k - self.target
                if context_start < 0 or context_end > len(self[chrom]):
                    yield None
                    continue
                else:
                    context = self[chrom][context_start:context_end].seq
            else:
                context_start = pos - self.k + self.target + 1
                context_end = pos + self.target + 1
                if context_start < 0 or context_end > len(self[chrom]):
                    yield None
                    continue
                else:
                    context = reverse_complement(self[chrom][context_start:
                                                             context_end].seq)
            if not re.match('^[ACGT]+$', context):
                context = None
            yield context

    def targets(self,
                bed: Union[str, TextIO] = None) -> Dict[str, int]:
        r"""return a dictionary of the number of sites of each k-mer

        Args:
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
