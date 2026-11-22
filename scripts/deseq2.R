suppressPackageStartupMessages({
    library(tximport)
    library(DESeq2)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Usage: deseq2.R <samplesheet.csv> <tx2gene.tsv> <output_dir>")
}

samplesheet_file <- args[1]
tx2gene_file <- args[2]
output_dir <- args[3]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)

# Read and validate metadata
samples <- read.csv(samplesheet_file, stringsAsFactors = FALSE)

required_columns <- c(
    "sample", "run_accession", "population",
    "temperature", "replicate"
)

if (!all(required_columns %in% colnames(samples))) {
    stop("Samplesheet is missing required columns")
}

if (anyDuplicated(samples$sample) || anyDuplicated(samples$run_accession)) {
    stop("Duplicate sample names or run accessions found")
}

samples$temperature <- factor(
    samples$temperature,
    levels = c(21, 29)
)

rownames(samples) <- samples$sample

if (any(table(samples$temperature) != 3)) {
    stop("Expected exactly three biological replicates per temperature")
}

# Locate Salmon quant.sf files by run accession
quant_files <- list.files(
    ".",
    pattern = "^quant\\.sf$",
    recursive = TRUE,
    full.names = TRUE
)

quant_accessions <- basename(dirname(quant_files))

file_index <- match(samples$run_accession, quant_accessions)

if (anyNA(file_index)) {
    stop(
        "Missing Salmon quantification for: ",
        paste(samples$run_accession[is.na(file_index)], collapse = ", ")
    )
}

files <- quant_files[file_index]
names(files) <- samples$sample

# Import transcript estimates and aggregate to genes
tx2gene <- read.delim(
    tx2gene_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

stopifnot(all(c("TXNAME", "GENEID") %in% colnames(tx2gene)))

txi <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene[, c("TXNAME", "GENEID")],
    countsFromAbundance = "no"
)

# Construct DESeq2 dataset
dds <- DESeqDataSetFromTximport(
    txi,
    colData = samples,
    design = ~ temperature
)

# Filter genes with very low counts
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Differential expression: 29C versus 21C
dds <- DESeq(dds)

res <- results(
    dds,
    contrast = c("temperature", "29", "21"),
    alpha = 0.05
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

# Export results
write.csv(
    res_df,
    file.path(output_dir, "differential_expression.csv"),
    row.names = FALSE
)

write.csv(
    counts(dds, normalized = FALSE),
    file.path(output_dir, "gene_counts_raw.csv")
)

write.csv(
    counts(dds, normalized = TRUE),
    file.path(output_dir, "gene_counts_normalized.csv")
)

write.csv(
    as.data.frame(colData(dds)),
    file.path(output_dir, "sample_metadata.csv")
)

# PCA plot
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(
    vsd,
    intgroup = "temperature",
    returnData = TRUE
)

percent_variance <- round(
    100 * attr(pca_data, "percentVar")
)

pca_plot <- ggplot(
    pca_data,
    aes(
        x = PC1,
        y = PC2,
        color = temperature,
        label = name
    )
) +
    geom_point(size = 4) +
    geom_text(vjust = -1, size = 3) +
    xlab(paste0("PC1: ", percent_variance[1], "% variance")) +
    ylab(paste0("PC2: ", percent_variance[2], "% variance")) +
    theme_bw() +
    labs(
        title = "RNA-seq PCA",
        color = "Temperature (°C)"
    )

ggsave(
    file.path(output_dir, "plots", "pca.png"),
    pca_plot,
    width = 8,
    height = 6,
    dpi = 300
)

# Volcano plot
volcano_data <- res_df
volcano_data$significance <- "Not significant"

volcano_data$significance[
    !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange > 0
] <- "Higher at 29°C"

volcano_data$significance[
    !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange < 0
] <- "Lower at 29°C"

volcano_plot <- ggplot(
    volcano_data,
    aes(
        x = log2FoldChange,
        y = -log10(pmax(padj, 1e-300)),
        color = significance
    )
) +
    geom_point(alpha = 0.65, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() +
    labs(
        title = "Differential expression: 29°C versus 21°C",
        x = "Log2 fold change",
        y = "-Log10 adjusted p-value",
        color = "Expression"
    )

ggsave(
    file.path(output_dir, "plots", "volcano.png"),
    volcano_plot,
    width = 8,
    height = 6,
    dpi = 300
)

# Record analysis details
writeLines(
    capture.output(sessionInfo()),
    file.path(output_dir, "sessionInfo.txt")
)

writeLines(
    c(
        paste("Genes tested:", nrow(res_df)),
        paste("Genes with adjusted p-value < 0.05:",
              sum(res_df$padj < 0.05, na.rm = TRUE)),
        paste("Genes higher at 29C:",
              sum(res_df$padj < 0.05 &
                  res_df$log2FoldChange > 0, na.rm = TRUE)),
        paste("Genes lower at 29C:",
              sum(res_df$padj < 0.05 &
                  res_df$log2FoldChange < 0, na.rm = TRUE))
    ),
    file.path(output_dir, "summary.txt")
)

message("DESeq2 analysis completed successfully.")
