include: "mutyper.smk"


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
