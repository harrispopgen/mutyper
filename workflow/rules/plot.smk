include: "mutyper.smk"


from itertools import combinations


rule plot_spectra:
    input:
        stratify=expand(rules.mutyper_spectra_stratify.output.spectra, rgn=RGNS),
        full=rules.mutyper_spectra.output.spectra,
    output:
        violin="results/plots/spectra/violin.pdf",
        heatmap="results/plots/spectra/heatmap.pdf",
        pca="results/plots/spectra/pca.pdf",
    log:
        "logs/plots/spectra.log",
    conda:
        "../envs/R.yml"
    script:
        "../scripts/spectra.R"


rule plot_spectra_stratify:
    input:
        full=rules.mutyper_spectra_stratify.output.spectra,
    output:
        violin="results/plots/spectra/stratify/{rgn}/violin.pdf",
        heatmap="results/plots/spectra/stratify/{rgn}/heatmap.pdf",
        pca="results/plots/spectra/stratify/{rgn}/pca.pdf",
    log:
        "logs/plots/spectra_{rgn}.log",
    conda:
        "../envs/R.yml"
    script:
        "../scripts/spectra.R"


rule plot_comparison:
    input:
        lambda wc: (rules.mutyper_spectra_stratify.output.spectra).format(rgn=wc.name1),
        lambda wc: (rules.mutyper_spectra_stratify.output.spectra).format(rgn=wc.name2),
    output:
        plot="results/plots/spectra/stratify/compare/all/{name1}_{name2}.pdf",
        fold="results/plots/spectra/stratify/compare/logfold/{name1}_{name2}.pdf",
    log:
        "logs/plots/compare/{name1}_{name2}.log",
    conda:
        "../envs/R.yml"
    script:
        "../scripts/compare.R"


def make_plot_comparison_outputs(wc):
    for name1, name2 in combinations(RGNS, 2):
        rtn = expand(rules.plot_comparison.output, name1=name1, name2=name2)
        for f in rtn:
            yield f


rule make_comparisons:
    input:
        make_plot_comparison_outputs,
