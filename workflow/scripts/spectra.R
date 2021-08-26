library(vroom)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(scales)
library(cowplot)
library(ggrepel)
library(data.table)
library(tidyverse)

spectra_f <- "https://eichlerlab.gs.washington.edu/help/mvollger/share/mutyper/spectra.txt"
spectra_f <- snakemake@input$full

heatmap_f <- "results/plots/spectra/heatmap.pdf"
violin_f <- "results/plots/spectra/heatmap.pdf"
pca_f <- "results/plots/spectra/heatmap.pdf"
heatmap_f <- snakemake@output$heatmap
pca_f <- snakemake@output$pca
violin_f <- snakemake@output$violin



make_spectra_matrix <- function(spectra) {
    spec.m <- as.matrix(spectra[, -1])
    row.names(spec.m) <- t(spectra[, 1])
    spec.m
}

make_spectra_long <- function(spectra) {
    long <- colnames(spectra)[grepl(".*>.*", colnames(spectra))]
    spectra %>%
        pivot_longer(
            cols = long,
            names_to = "spectra",
            values_to = "count"
        ) %>%
        mutate(first_two_bases = substr(spectra, 1, 2)) %>%
        data.table()
}


read_spectra <- function(f) {
    spectra <- fread(spectra_f, sep = "\t")
    spectra.m <- make_spectra_matrix(spectra)
    pca_res <- prcomp(spectra.m, center = TRUE, scale. = TRUE)
    spectra$PC1 <- pca_res$x[, 1]
    spectra$PC2 <- pca_res$x[, 2]

    list(
        df = spectra,
        long = make_spectra_long(spectra),
        pca = pca_res,
        m = spectra.m
    )
}
spectra <- read_spectra(spectra_f)

pdf(heatmap_f, height = 8, width = 8)
heatmap(spectra$m)
dev.off()

pdf(pca_f, height = 12, width = 12)
autoplot(spectra$pca, data = spectra$df) +
    geom_text_repel(aes(label = sample)) +
    theme_cowplot() +
    theme(legend.position = "top")
dev.off()


p.violin <- spectra$long %>%
    ggplot(aes(x = spectra, y = count, color = first_two_bases)) +
    geom_violin() +
    geom_jitter() +
    scale_x_discrete(guide = guide_axis(n.dodge = 6)) +
    cowplot::theme_minimal_vgrid() +
    theme(legend.position = "none")

scale <- 1
ggsave(violin_f, height = 12 * scale, width = 24 * scale, plot = p.violin)