library(vroom)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(scales)
library(cowplot)
library(ggrepel)
library(data.table)
library(tidyverse)
library(ggforce)
library(glue)

spectra_f_1 <- "results/spectra/stratify/SD_spectra.txt"
spectra_f_2 <- "results/spectra/stratify/Unique_spectra.txt"
spectra_f_1 <- snakemake@input[[1]]
spectra_f_2 <- snakemake@input[[2]]

name1 <- "SD"
name2 <- "Unique"
name1 <- snakemake@wildcards$name1
name2 <- snakemake@wildcards$name2

out_plot <- "~/Desktop/1_2.pdf"
out_fold <- "~/Desktop/log_fold.pdf"
out_fold <- snakemake@output$fold
out_plot <- snakemake@output$plot



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
    spectra <- fread(f, sep = "\t")
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
spec1 <- read_spectra(spectra_f_1)
spec2 <- read_spectra(spectra_f_2)
l <- list(spec1$long, spec2$long)
names(l) <- c(name1, name2)
spec <- bind_rows(l, .id = "stratify") %>%
    group_by(stratify) %>%
    mutate(percent = 100 * count / sum(count)) %>%
    data.table()

pdf(out_plot, height = 11, width = 8)
for (two in unique(spec$first_two_bases)) {
    p <- spec %>%
        filter(first_two_bases == two) %>%
        ggplot(
            aes(
                x = spectra,
                y = percent,
                color = spectra,
            )
        ) +
        geom_violin() +
        geom_jitter() +
        scale_x_discrete(guide = guide_axis(n.dodge = 6)) +
        cowplot::theme_minimal_vgrid() +
        facet_grid(~stratify) +
        theme(legend.position = "none")
    print(p)
}
dev.off()


fold_change.df <- spec %>%
    group_by(stratify, spectra) %>%
    summarise(count = sum(count)) %>%
    ungroup() %>%
    group_by(stratify) %>%
    mutate(percent = 100 * count / sum(count)) %>%
    ungroup() %>%
    pivot_wider(
        id_cols = spectra,
        names_from = stratify,
        values_from = "percent"
    ) %>%
    mutate(
        log_change = log2((!!as.name(name1)) / (!!as.name(name2)))
    ) %>%
    drop_na() %>%
    data.table()


p.fold <- fold_change.df %>%
    ggplot(
        aes(
            x = spectra,
            y = log_change,
            color = log_change > 0,
            fill = log_change > 0
        )
    ) +
    ylab(glue("log2({name1}/{name2})")) +
    geom_bar(stat = "identity") +
    scale_x_discrete(guide = guide_axis(n.dodge = 6)) +
    cowplot::theme_minimal_vgrid() +
    theme(legend.position = "none")
scale <- 1
ggsave(out_fold, height = 12 * scale, width = 24 * scale, plot = p.fold)