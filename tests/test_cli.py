#! /usr/bin/env python

import pytest
from mutyper import cli
import argparse
import pandas as pd
import io


def test_spectra(capsys):
    args = argparse.Namespace(
        vcf="tests/test_data/snps.vcf", population=False, randomize=False
    )
    cli.spectra(args)
    captured = capsys.readouterr()
    df = pd.read_csv(io.StringIO(captured.out), sep="\t", index_col=0)
    df_target = pd.DataFrame(
        {"ACA>ATA": [0, 1], "ACC>ATC": [1, 0], "ACG>ATG": [1, 1], "ACT>ATT": [2, 0]},
        index=pd.Index(["sample1", "sample2"], name="sample"),
    )
    pd.testing.assert_frame_equal(df, df_target)


def test_spectra_randomize(capsys):
    args = argparse.Namespace(
        vcf="tests/test_data/snps.vcf", population=False, randomize=True
    )
    cli.spectra(args)
    captured = capsys.readouterr()
    df = pd.read_csv(io.StringIO(captured.out), sep="\t", index_col=0)
    # there are two possibilities due to haplotype randomization
    df_target1 = pd.DataFrame(
        {"ACA>ATA": [0, 1], "ACC>ATC": [1, 0], "ACG>ATG": [1, 0], "ACT>ATT": [1, 0]},
        index=pd.Index(["sample1", "sample2"], name="sample"),
    )
    df_target2 = pd.DataFrame(
        {"ACA>ATA": [0, 1], "ACC>ATC": [1, 0], "ACG>ATG": [0, 1], "ACT>ATT": [1, 0]},
        index=pd.Index(["sample1", "sample2"], name="sample"),
    )

    try:
        pd.testing.assert_frame_equal(df, df_target1)
    except AssertionError:
        pd.testing.assert_frame_equal(df, df_target2)


def test_spectra_haploid(capsys):
    args = argparse.Namespace(
        vcf="tests/test_data/snps.haploid.vcf", population=False, randomize=False
    )
    cli.spectra(args)
    captured = capsys.readouterr()
    df = pd.read_csv(io.StringIO(captured.out), sep="\t", index_col=0)
    df_target = pd.DataFrame(
        {"ACA>ATA": [0, 1], "ACC>ATC": [0, 0], "ACG>ATG": [1, 1], "ACT>ATT": [1, 0]},
        index=pd.Index(["sample1", "sample2"], name="sample"),
    )
    pd.testing.assert_frame_equal(df, df_target)


def test_spectra_missing_gts(capsys, caplog):
    args = argparse.Namespace(
        vcf="tests/test_data/snps.missing_gts.vcf", population=False, randomize=False
    )
    cli.spectra(args)
    captured = capsys.readouterr()
    df = pd.read_csv(io.StringIO(captured.out), sep="\t", index_col=0)
    df_target = pd.DataFrame(
        {"ACA>ATA": [0, 1], "ACC>ATC": [1, 0], "ACG>ATG": [1, 1], "ACT>ATT": [1, 1]},
        index=pd.Index(["sample1", "sample2"], name="sample"),
    )
    pd.testing.assert_frame_equal(df, df_target)
    assert "Ambiguous genotypes found" in caplog.text


def test_variant_missing_gts_nonstrict(capsys):
    args = argparse.Namespace(
        fasta="tests/test_data/ancestor.fa", vcf="tests/test_data/snps.missing_gts.variants.vcf", k=3, target=None, sep=":", chrom_pos=0, strand_file=None, strict=False
    )
    cli.variants(args)
    captured = capsys.readouterr()
    df = pd.read_csv(io.StringIO(captured.out), skiprows=8, sep="\t")
    df_info_expand = df["INFO"].str.split(";", expand=True)
    df_info_expand.columns = ["AN", "AC", "AF", "mutation_type"]
    df_subset = pd.concat([df[["POS", "REF", "ALT", "sample1", "sample2"]], df_info_expand], axis=1)
    df_target = pd.DataFrame(
        {
            "POS": [2,3,4,7,8,11],
            "REF": ["A", "A", "C", "C", "C", "A"],
            "ALT": ["T", "C", "T", "T", "T", "G"],
            "sample1": ["1|.", ".|0", ".|.", "0|1", "1|1", "0/0"],
            "sample2": ["1/0", "1/1", "0/1", "1/0", "1|1", "0/0"],
            "AN": [f"AN={an}" for an in [3, 3, 2, 4, 4, 4]],
            "AC": [f"AC={ac}" for ac in [2, 2, 1, 2, 4, 0]],
            "AF": [f"AF={af}" for af in [0.67, 0.666667, 0.5, 0.50, 1, 0]],
            "mutation_type": [f"mutation_type={type}" for type in ["AAA>ATA", "AAC>ACC", "ACC>ATC", "CCG>CTG", "CCC>CTC", "AAA>AGA"]]
        }
    )
    pd.testing.assert_frame_equal(df_subset, df_target)


def test_variant_missing_gts_strict(capsys):
    args = argparse.Namespace(
        fasta="tests/test_data/ancestor.fa", vcf="tests/test_data/snps.missing_gts.variants.vcf", k=3, target=None, sep=":", chrom_pos=0, strand_file=None, strict=True
    )
    cli.variants(args)
    captured = capsys.readouterr()
    df = pd.read_csv(io.StringIO(captured.out), skiprows=8, sep="\t")
    df_info_expand = df["INFO"].str.split(";", expand=True)
    df_info_expand.columns = ["AN", "AC", "AF", "mutation_type"]
    df_subset = pd.concat([df[["POS", "REF", "ALT", "sample1", "sample2"]], df_info_expand], axis=1)
    df_strict_target = pd.DataFrame(
        {
            "POS": [2,3,4,11],
            "REF": ["A", "A", "C", "A"],
            "ALT": ["T", "C", "T", "G"],
            "sample1": ["1|.", ".|0", ".|.", "0/0"],
            "sample2": ["1/0", "1/1", "0/1", "0/0"],
            "AN": [f"AN={an}" for an in [3, 3, 2, 4]],
            "AC": [f"AC={ac}" for ac in [2, 2, 1, 0]],
            "AF": [f"AF={af}" for af in [0.67, 0.666667, 0.5, 0]],
            "mutation_type": [f"mutation_type={type}" for type in ["AAA>ATA", "AAC>ACC", "ACC>ATC", "AAA>AGA"]]
        }
    )
    pd.testing.assert_frame_equal(df_subset, df_strict_target)


def test_ksfs(capsys):
    args = argparse.Namespace(vcf="tests/test_data/snps.vcf", k=3)
    cli.ksfs(args)
    captured = capsys.readouterr()
    df = pd.read_csv(io.StringIO(captured.out), sep="\t", index_col=0)
    df_target = pd.DataFrame(
        {
            "ACA>ATA": [1, 0, 0],
            "ACC>ATC": [1, 0, 0],
            "ACG>ATG": [0, 1, 0],
            "ACT>ATT": [0, 1, 0],
        },
        index=pd.Index([1, 2, 3], name="sample_frequency"),
    )
    pd.testing.assert_frame_equal(df, df_target)


def test_ksfs_missing_gts():
    args = argparse.Namespace(vcf="tests/test_data/snps.missing_gts.vcf", k=3)
    with pytest.raises(
        ValueError, match=r"different AN [0-9]* and [0-9]* indicates missing genotypes"
    ):
        cli.ksfs(args)
