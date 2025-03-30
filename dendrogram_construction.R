library(rafalib) #nice plot arrangement

# Load the airway dataset
# Assuming `airway` is the SummarizedExperiment object
counts <- assay(airway)  # Extract raw counts
metadata <- colData(airway) %>% as.data.frame()  # Extract sample metadata

#Convert it to a distance object
d <- dist( t(counts) , method="euclidean")

#Compute sample correlations
sample_cor <- cor( counts )
round(sample_cor,4)
pheatmap::pheatmap(sample_cor)

#Transform the scale from correlations
cor_distance <- -(sample_cor-1)/2
round(cor_distance,4)
pheatmap::pheatmap(cor_distance)

#Clustering using euclidean distance
{
  mypar(1,2,mar=c(6,4,2,1))
  h <- hclust(d,method="complete")
  plot( as.dendrogram(h),las=1,main="d=euclidean\nh=complete")
  points(1:ncol(counts),rep(0,ncol(counts)),pch=16,cex=2,col=metadata$dex[h$order])
}


h2 <- hclust(d,method="complete")
{
  plot( as.dendrogram(h2),las=1, main="d=correlation\nh=complete")
  points(1:ncol(counts),rep(0,ncol(counts)),pch=16,cex=2, col=metadata$dex[h2$order])
}

h1 <- hclust(d,method="complete")
h2 <- hclust(d,method="complete")

# Convert to dendrogram
dend <- as.dendrogram(h2)
# Color the dendrogram labels based on metadata
labels_colors(dend) <- metadata$dex[h2$order]
# Convert to a ggdendro-compatible format
dend_data <- ggdendro::dendro_data(dend)

# Convert to dendrogram
dend <- as.dendrogram(h1)
# Color the dendrogram labels based on metadata
labels_colors(dend) <- metadata$dex[h1$order]
# Convert to a ggdendro-compatible format
dend_data1 <- ggdendro::dendro_data(dend)

# Create ggplot-based dendrogram
ggplot() +
  geom_segment(data = dend_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = data.frame(x = 1:ncol(counts), y = 0, color = metadata$dex[h2$order]),
             aes(x = x, y = y, color = color), size = 4) +
  scale_color_manual(values = unique(metadata$dex), name = "Metadata Group") +
  theme_minimal() +
  labs(title = "Complete") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggplot() +
  geom_segment(data = dend_data1$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = data.frame(x = 1:ncol(counts), y = 0, color = metadata$dex[h1$order]),
             aes(x = x, y = y, color = color), size = 4) +
  scale_color_manual(values = unique(metadata$dex), name = "Metadata Group") +
  theme_minimal() +
  labs(title = "Correlation") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
