# Load necessary libraries
# Assume 'vst_counts' is a matrix or data frame containing VST-transformed RNA-seq counts
# Rows are genes, columns are samples

# Step 1: Select the 500 most variable genes
top_variable_genes <- head(order(apply(assay(vstcounts), 1, var), decreasing = TRUE), 500) #assay(vstcounts) or assay(vst(dds, blind = T)) or normalized cts

# Subset the data to keep only these genes
vst_top500 <- assay(vstcounts)[top_variable_genes, ]

# Step 2: Compute Spearman correlation distance
spearman_dist <- as.dist(1 - cor(vst_top500, method = "spearman"))

# Perform hierarchical clustering using ward.D2 method
sample_clustering <- hclust(spearman_dist, method = "ward.D2")

# Plot the dendrogram
plot(sample_clustering, main = "Hierarchical Clustering of Samples", xlab = "", sub = "")

# Optional: Generate a heatmap with hierarchical clustering
pheatmap(vst_top500, clustering_distance_cols = spearman_dist, clustering_method = "ward.D2")
