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
enrichment_matrix_test <- performGOEnrichment( genes_check, go_terms_check, p.adjust_method = 'BH')
test_hclust_results <- performHierarchicalClustering( enrichment_matrix_test )
test_clusters <- groupClusters( test_hclust_results)
test_cluster_labels <- annotateClusters( test_clusters )

go_terms_check[["TEST"]] <- unique(rapply(go_terms_check, function(x) head(x,Inf)))


# Heatmaps creation
svg( "test_heatmap.svg", width = 16, height = 8 )
plotHeatmap( enrichment_matrix_test, test_hclust_results, test_clusters, test_cluster_labels, min_num_terms = 1, heatmap_colors = 'auto',
             dendrogram_lwd = 1, header_lwd = 1, cluster_label_cex = 1, sample_label_cex = 1, cluster_label_width = 0.8)
dev.off()

png( "test_heatmap.png", width = 9600, height = 7200 )
plotHeatmap( enrichment_matrix_test, test_hclust_results, test_clusters, test_cluster_labels, min_num_terms = 1, heatmap_colors = 'auto',
             dendrogram_lwd = 8, header_lwd = 8, cluster_label_cex = 8, sample_label_cex = 8, cluster_label_width = 0.8)
dev.off()


# goSTAG analysis
source("https://bioconductor.org/biocLite.R")
biocLite("goSTAG")
library(goSTAG)
library(biomaRt)
gos <- loadGeneLists( "test/test_data/example_go_terms.gmt" )
genes <- loadGeneLists("test/test_data/example_DE_gene.gmt" )

# Add all gene ID in "ALL"
query <- read.csv(file="test/test_data/example_gene_go_query.tsv", header=TRUE, sep="\t")
query = query[ query[,1]!="", ]
query = query[ query[,2]!="", ]

## Filter the BioMart query to only include valid GO terms from the specified domain
go_ids = unique(query[,2])

domains = Ontology(GOTERM)

go_ids = unique(query[,2])
go_ids_filtered = go_ids[domains[names(go_ids)]==toupper('BP') ]
go_ids_filtered = go_ids_filtered[ ! is.na(go_ids_filtered) ]
query_filtered = query[ query[,2] %in% go_ids_filtered, ]

## Get the list of gene annotations for each GO term
go_genes = lapply( go_ids_filtered, function(x) unique(query_filtered[query_filtered[,2]==x,1]) )
names(go_genes) = go_ids_filtered

## Add full list of all annotated genes
gos[[ "ALL" ]] = unique(query[,1])

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
