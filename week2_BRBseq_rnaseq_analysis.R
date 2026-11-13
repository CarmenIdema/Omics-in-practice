library(DESeq2)
library(data.table)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

counts_df <- fread("data/gene_counts.tsv", data.table = FALSE, check.names = FALSE)
metadata <- fread("metadata/sample_metadata.tsv", data.table = FALSE)
gene_ids <- counts_df[[1]]
counts <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(counts) <- gene_ids
mode(counts) <- "numeric"

stopifnot(!anyDuplicated(rownames(counts)))
stopifnot(!anyDuplicated(metadata$sample))
stopifnot(all(counts >= 0))
stopifnot(all(counts == round(counts)))
stopifnot(setequal(colnames(counts), metadata$sample))
metadata <- metadata[match(colnames(counts), metadata$sample), ]
stopifnot(identical(colnames(counts), metadata$sample))
rownames(metadata) <- metadata$sample

dim(counts)
summary(colSums(counts))
table(metadata$condition)
head(metadata)

lib_df <- data.frame(sample = colnames(counts),
library_size = colSums(counts),
condition = metadata$condition)
p_lib <- ggplot(lib_df, aes(sample, library_size, fill = condition)) +
geom_col() +
scale_y_continuous(labels = scales::comma) +
labs(x = NULL, y = "Total raw counts", title = "Library sizes") +
theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/library_sizes.pdf", p_lib, width = 8, height = 5)

metadata$condition <- factor(metadata$condition)
metadata$condition <- relevel(metadata$condition, ref = "K562_0h")

dds <- DESeqDataSetFromMatrix(
countData = counts,
colData = metadata,
design = ~ condition
)

# Example only; confirm with the instructor
# metadata$batch <- factor(metadata$batch)
# dds <- DESeqDataSetFromMatrix(counts, metadata, design = ~ batch + condition)

keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
cat("Genes retained:", nrow(dds), "\n")

saveRDS(dds, "results/dds_prefiltered.rds")

vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

pca_data <- plotPCA(vsd, intgroup = c("condition", "replicate"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, shape = replicate, label = name)) +
geom_point(size = 3) +
geom_text_repel(show.legend = FALSE, max.overlaps = Inf) +
xlab(paste0("PC1: ", percent_var[1], "% variance")) +
ylab(paste0("PC2: ", percent_var[2], "% variance")) +
coord_equal() + theme_bw() +
labs(title = "PCA of variance-stabilized counts")
ggsave("figures/PCA_samples.pdf", p_pca, width = 7, height = 5.5)
fwrite(pca_data, "results/PCA_coordinates.tsv", sep = "\t")

sample_dists <- dist(t(vsd_mat))
sample_dist_mat <- as.matrix(sample_dists)
rownames(sample_dist_mat) <- colnames(vsd_mat)
colnames(sample_dist_mat) <- colnames(vsd_mat)
ann_col <- metadata[, c("condition", "replicate"), drop = FALSE]
pdf("figures/sample_distance_heatmap.pdf", width = 7, height = 6)
pheatmap(sample_dist_mat, annotation_col = ann_col, annotation_row = ann_col,
main = "Sample-to-sample distances")
dev.off()
