#!/usr/bin Rscript

"
This script makes a goSTAG analysis using the go_gene.gmt and the DE_gene.gmt files.
go_gene.gmt contains the GO terms and all the genes they annotate.
DE_gene.gmt contains the DE gene lists and all the genes they contain.
These files are create by gmt_DE_gene_creation.py and gmt_go_genome_creation.py
"

args <- commandArgs(TRUE)
gene_gmt <- args[1]
go_gmt <- args[2]

# Test data.
go_terms_check <- loadGOTerms()
data( goSTAG_example_gene_lists )
genes_check <- head( lapply( goSTAG_example_gene_lists, head ) )

# goSTAG analysis
source("https://bioconductor.org/biocLite.R")
biocLite("goSTAG")
library(goSTAG)
gos <- loadGeneLists( "test/test_data/example_go_terms.gmt" )
genes <- loadGeneLists("test/test_data/example_DE_gene.gmt" )
gos[["ALL"]] <- unique(rapply(gos, function(x) head(x,Inf)))
enrichment_matrix <- performGOEnrichment( genes, gos, p.adjust_method = 'BH')
hclust_results <- performHierarchicalClustering( enrichment_matrix )
clusters <- groupClusters( hclust_results)
cluster_labels <- annotateClusters( clusters )

# Heatmaps creation
svg( "heatmap.svg", width = 16, height = 8 )
plotHeatmap( enrichment_matrix, hclust_results, clusters, cluster_labels, min_num_terms = 1, heatmap_colors = 'auto',
             dendrogram_lwd = 1, header_lwd = 1, cluster_label_cex = 1, sample_label_cex = 1, cluster_label_width = 0.8)
dev.off()

png( "heatmap.png", width = 9600, height = 7200 )
plotHeatmap( enrichment_matrix, hclust_results, clusters, cluster_labels, min_num_terms = 1, heatmap_colors = 'auto',
             dendrogram_lwd = 8, header_lwd = 8, cluster_label_cex = 8, sample_label_cex = 8, cluster_label_width = 0.8)
dev.off()
