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
