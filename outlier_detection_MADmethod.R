#using shinyngs function and the same has been implemented in nfcore tool
devtools::install_github('pinin4fjords/shinyngs', upgrade_dependencies = FALSE)
library(shinyngs)

plotdata = madScore(matrix = normalized_cts, #assay(vstcounts) or assay(vst(dds, blind = T)) or normalized cts
         sample_sheet = coldata %>% column_to_rownames("sample"), #rownames should match the column names of the matrix
         groupby = "tissue", #or condition of interest
         )

#scatterplot
ggplot(plotdata %>% rownames_to_column("sample"), aes(x = group, y = mad, color = outlier, label = sample)) +
  geom_point(size = 3, color = "orange") +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 10, color = "purple") +
  #geom_hline(yintercept = threshold, color = "orange", linetype = "dashed", size = 2) +
  labs(
    title = "MAD Scores by Sample Group",
    x = "Sample Group",
    y = "MAD Score",
    color = "Outlier"
  )
