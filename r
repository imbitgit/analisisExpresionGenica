if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")

if (!requireNamespace("pheatmap", quietly = TRUE))
    install.packages("pheatmap")

library(DESeq2)
library(ggplot2)
library(pheatmap)


set.seed(123)
counts <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
colData <- data.frame(
    condition = factor(rep(c("control", "treatment"), each = 5))
)
rownames(counts) <- paste0("gene", 1:100)
colnames(counts) <- paste0("sample", 1:10)

#Crear objeto DESeqDataSet y ejecutar análisis
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

#Graficar Volcano Plot
res$significant <- ifelse(res$padj < 0.05, "Significant", "Not Significant")

volcano <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("black", "red")) + 
    theme_minimal() +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")

print(volcano)

#Graficar Heatmap


if (!require("pheatmap")) {
  install.packages("pheatmap", repos = "http://cran.us.r-project.org")
  library(pheatmap)
}

#Normalizar los datos
vsd <- vst(dds, blind = FALSE)

top_genes <- head(order(res$padj, na.last = NA), 20)

if (length(top_genes) < 5) {
    message("Pocos genes significativos encontrados. Seleccionando por varianza para el gráfico.")
    #Calcular varianza de cada gen
    gene_vars <- apply(assay(vsd), 1, var)
    #Seleccionar los 20 con mayor varianza
    top_genes <- head(order(gene_vars, decreasing = TRUE), 20)
}

heatmap_data <- assay(vsd)[top_genes, ]

#Crear el mapa de calor
pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = TRUE, 
         show_colnames = TRUE,
         annotation_col = colData)
