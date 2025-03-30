#creating empty lists
pcaplots_list <- list()
screeplots_list <- list()
volcanoplots_list <- list()
#deg_results5 = list() #ALL BATS
#deg_results5 <- qs::qread("~/Bulk_RNAseq/NHCS_bulkRNAseq/04_data_objects/degresults5_ALLbats_vs_humanGTEx.qs")
#deg_results6 <- list() #Only PAR
deg_results_lv <- list() #deg_results_lv1 (Partial) or deg_results_lv (complete)
#deg_results_lv <- qs::qread("~/Bulk_RNAseq/NHCS_bulkRNAseq/04_data_objects/DEGresultslv_ALLcomparisons.qs") #Complete
loadings_plotlist <- list()
metadata_list <- list()
tissue_wise_anno <- list()
#deg_results_tfs1 = list() #for storing the results
#deg_results11 = list() #storing the deg results

comparisons2 = tidyr::crossing(group1 = unique(counts_completetpm_meta$species),
                group2 = unique(counts_completetpm_meta$species)) %>% 
  filter(row_number() %in% c(5,7,8,9)) #%>%
  #filter(row_number() %in% c(1,4)) #comment it out if not using ESPL data

for(i in 1:nrow(comparisons2)) { #[-c(3,4,7)]
  
  #cat("Comparing ALL and human transcriptomes for", i, "\n")
  cat(paste("Comparing", comparisons2 %>% filter(row_number() == i) %>% pull(group1), "vs", 
        comparisons2 %>% filter(row_number() == i) %>% pull(group2)), "\n")
  
  #subset the expression matrix to include only samples from that tissue
  bat_samples <- counts_completetpm_meta %>% 
                    filter(species == comparisons2 %>% filter(row_number() == i) %>% pull(group2)) %>% 
                    #filter(!samples %in% c("SRX24864131", "SRX24865093", "SRX9833474", "SRX9833464")) %>%
                    pull(samples)
  
  if(comparisons2 %>% filter(row_number() == i) %>% pull(group1) == "human") {
    
    hs_samples <- counts_completetpm_meta %>% # counts_completetpm_meta (complete) or gtex_meta (Partial)
        dplyr::filter(tissue == "heart" & SMTSD == "Heart - Left Ventricle") %>% # & SMTSD == "Heart - Left Ventricle" #Uncomment for the full dataset
        filter(species == comparisons2 %>% filter(row_number() == i) %>% pull(group1)) %>% #Uncomment for the full dataset
        #slice_max(SMRIN,n = (length(bat_samples)),with_ties = FALSE) %>% #Uncomment for the partial dataset
        pull(samples) #samples (complete) or SAMPID (subset)
    
    } else {
    
    hs_samples <- counts_completetpm_meta %>%
    filter(species == comparisons2 %>% filter(row_number() == i) %>% pull(group1)) %>%
    pull(samples)    
      
    }
  
  meta6 <- counts_completetpm_meta %>% filter(samples %in% c(hs_samples, bat_samples))
  meta6$species <- forcats::fct_relevel(meta6$species, comparisons2 %>% filter(row_number() == i) %>% pull(group2))
  
  metadata_list[[i]] <- data.frame(tissue = "heart",
                                   num_samples = length(hs_samples)
                                   )
  
  merged_tpm <- counts_completetpm %>% 
    dplyr::select(meta6$samples) #%>% 
    #dplyr::select(any_of(meta6$samples)) #%>% #same order as metadata
    #filter(rownames(.) %in% gene_filter) #remove it if no gene filter required
    
  print(all(colnames(merged_tpm) == rownames(meta6 %>% column_to_rownames("samples"))))
  
  #convert NA to 0 in the merged TPM matrix
  merged_tpm[is.na(merged_tpm)] = 0
  merged_tpm <- merged_tpm %>% mutate_if(is.numeric, list(round))
  
  ##PCA on log transform TPM or TPM counts using Neha' method- removed it from the code
  #trying PCA plot using PCAtools
  p <- PCAtools::pca(log2(merged_tpm + 1), metadata = meta6 %>% column_to_rownames("samples"), removeVar = 0.1)
  p2 <- biplot(p, lab = paste0(p$metadata$species), colby = 'species',hline = 0, vline = 0,legendPosition = 'right')
  #pcaplots_list[[i]] <- p2
  #biplot(p, showLoadings = TRUE, lab = NULL)
  #plotloadings(p,rangeRetain = 0.01,labSize = 4.0,title = 'Loadings plot',subtitle = 'PC1, PC2, PC3, PC4, PC5',caption = 'Top 1% variables',shape = 24,col = c('limegreen', 'black', 'red3'),drawConnectors = TRUE)
  #p1 <- plotloadings(p, components = getComponents(p, 1:5),rangeRetain = 0.1,labSize = 4.0,absolute = FALSE,title = paste('Loadings plot for', i), subtitle = 'PCs',caption = 'Top 10% variables',shape = 23, shapeSizeRange = c(1, 16),col = c('white', 'pink'),drawConnectors = FALSE)
  #loadings_plotlist[[i]] <- p1
  q <- screeplot(p, axisLabSize = 18, titleLabSize = 22, title = paste("SCREE plot for", i))
  screeplots_list[[i]] <- q 
  #biplot(p, showLoadings = TRUE, labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
  
  
  #DGE analysis using Wilcoxon test
  pvalues <- sapply(1:nrow(merged_tpm),function(j){
    data <- cbind.data.frame(gene=as.numeric(t(merged_tpm[j,])),species = factor(meta6$species))
    p = wilcox.test(gene ~ species, data)$p.value
    return(p)
  })
     
  fdr=p.adjust(pvalues,method = "fdr")   
  
  #### Calculate fold-change for each gene
  conditionsLevel <- levels(factor(meta6$species))
  dataCon1=merged_tpm[,c(which(factor(meta6$species)==conditionsLevel[1]))]
  dataCon2=merged_tpm[,c(which(factor(meta6$species)==conditionsLevel[2]))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  
  #### Output results base on FDR threshold
  outRst <- data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
  rownames(outRst) = rownames(merged_tpm)
  outRst = na.omit(outRst)
  
  print(outRst %>% filter_all(all_vars(!is.infinite(.))) %>% 
    dplyr::mutate(DEG = case_when(log2foldChange > 0.5 & FDR < 0.01 ~ paste0(comparisons2 %>% filter(row_number() == i) %>% pull(group1),"-enriched"), #enriched in Humans
                               log2foldChange < -0.5 & FDR < 0.01 ~ paste0(comparisons2 %>% filter(row_number() == i) %>% pull(group2),"-enriched"), #enriched in Bats
                               TRUE ~ "non-sig")) %>% 
  dplyr::count(DEG))
  
  deg_results_lv[[paste0(comparisons2 %>% filter(row_number() == i) %>% pull(group1),"_vs_",comparisons2 %>% filter(row_number() == i) %>% pull(group2))]] <- outRst %>% filter_all(all_vars(!is.infinite(.))) %>% 
    dplyr::mutate(DEG = case_when(log2foldChange > 0.5 & FDR < 0.01 ~ "up", #enriched in Humans
                                  log2foldChange < -0.5 & FDR < 0.01 ~ "down", #enriched in Bats
                                  TRUE ~ "non-sig"),
                  enrichment = case_when(DEG == "up" ~ paste0(comparisons2 %>% filter(row_number() == i) %>% pull(group1),"-enriched"),
                                         DEG == "down" ~ paste0(comparisons2 %>% filter(row_number() == i) %>% pull(group2),"-enriched"),
                                         TRUE ~ "notsig"
                                         ),
                  tissue = "heart",
                  comparison = paste0(comparisons2 %>% filter(row_number() == i) %>% pull(group1),"_vs_",comparisons2 %>% filter(row_number() == i) %>% pull(group2)))
  
  r <- EnhancedVolcano::EnhancedVolcano(outRst %>% filter_all(all_vars(!is.infinite(.))),
                  lab = outRst %>% filter_all(all_vars(!is.infinite(.))) %>% rownames(),
                  x = 'log2foldChange',
                  y = 'FDR',
                  pCutoff = 0.01,
                  FCcutoff = 1,                
                  title = paste("Comparison of bulk transcriptomes between human and bat", i),
                  caption = "|Log2FC| > 1; padj < 0.01",
                  subtitle = "Differential expression analysis",
                  #pointSize = 3.0,
                  #colAlpha = 0.7,
                  #legendLabels=c('Not sig.','Log2FC','padj',
                   #              'padj & Log2FC'),
                  #drawConnectors = TRUE,
                  labSize = 4.0
                  )
  
  volcanoplots_list[[i]] <- r 
  
}

#create patchworks
#+ plot_annotation(title = 'Volcano plots', subtitle = 'human-enriched DEGs (log2FC > 0) & bat-enriched DEGs (log2FC < 0)')
pcaplots_list[1]
loadings_plotlist[1]
bind_rows(metadata_list)

qs::qsave(deg_results_lv, "~/Bulk_RNAseq/NHCS_bulkRNAseq/04_data_objects/DEGresultslv_ALLcomparisons.qs")

#save excel files
#writexl::write_xlsx(list_rbind(map(deg_results_lv1, function(x) {x %>% rownames_to_column("gene")}), names_to = "comparison"), "~/Bulk_RNAseq/NHCS_bulkRNAseq/07_csv_outputs/DEGresults_Balanced_LVGTEx_vs_others_heartcomparisons.xlsx")
