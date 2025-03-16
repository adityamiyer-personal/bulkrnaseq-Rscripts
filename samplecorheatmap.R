### Sample correlation heatmap using vst counts
sampleCor <- cor(normalized_cts)
sampleDist <- as.dist(1 - cor(normalized_cts))
sampleDistMatrix <- as.matrix(sampleDist)
blueColours <- brewer.pal(9, "Blues")
colors <- colorRampPalette(rev(blueColours))(255)
annotation_row = data.frame(
                    Group = str_before_first(colnames(sampleDistMatrix), "_")
                )
rownames(annotation_row) <- rownames(sampleDistMatrix)
pheatmap(sampleDistMatrix,
                   clustering_distance_cols = sampleDist, 
                   color = colors,
                   fontsize_row = 6,
                   fontsize_col = 8,
                   annotation_row = annotation_row,
                   cutree_rows = 4
                   )
