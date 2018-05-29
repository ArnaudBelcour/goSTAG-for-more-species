#!/usr/bin Rscript

"
This script makes a goSTAG analysis using the go_gene.gmt and the DE_gene.gmt files.
go_gene.gmt contains the GO terms and all the genes they annotate.
DE_gene.gmt contains the DE gene lists and all the genes they contain.
These files are create by gmt_DE_gene_creation.py and gmt_go_genome_creation.py
"

# goSTAG analysis
source("https://bioconductor.org/biocLite.R")
biocLite("goSTAG")
library(goSTAG)
library(biomaRt)
library(GO.db)
gos <- loadGeneLists( "go_gene.gmt" )
genes <- loadGeneLists("DE_gene.gmt" )

# Add all gene ID in "ALL"
query <- read.csv(file="query_go.tsv", header=TRUE, sep="\t")
query = query[ query[,1]!="", ]
query = query[ query[,2]!="", ]

## Filter the BioMart query to only include valid GO terms from the specified domain
go_ids = unique(query[,2])

domains = Ontology(GOTERM)

go_ids = unique(query[,2])
go_ids_filtered = go_ids[domains[names(go_ids)]==toupper('MF') ]
go_ids_filtered = go_ids_filtered[ ! is.na(go_ids_filtered) ]
query_filtered = query[ query[,2] %in% go_ids_filtered, ]

## Get the list of gene annotations for each GO term
go_genes = lapply( go_ids_filtered, function(x) unique(query_filtered[query_filtered[,2]==x,1]) )
names(go_genes) = go_ids_filtered

## Add full list of all annotated genes
gos[[ "ALL" ]] = unique(query[,1])
gos[["ALL"]] <- unique(rapply(gos, function(x) head(x,Inf)))

enrichment_matrix <- performGOEnrichment( genes, gos)
hclust_results <- performHierarchicalClustering( enrichment_matrix, distance_method = "euclidean", clustering_method = "complete" )
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

