# pcaData <- plotPCA(rld_pk_filtered_coding, intgroup = c("kit_name", "filter_name"), returnData = TRUE, ntop = 500)


# calculate the variance for each gene
rv <- rowVars(assay(rld_pk_filtered_coding))

# Top n genes by variance to keep.
ntop <- 500

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- stats::prcomp(t(assay(rld_pk_filtered_coding)[select,]))

# Loadings for the first two PCs.
loadings <- pca$rotation[, seq_len(2)]

loadings_df <- as.data.frame(loadings) %>% 
  # mutate(gene_id = rownames(loadings))
  tibble::rownames_to_column(var = "gene_id") %>% 
  left_join(gene_annotation, by = "gene_id")

loadings_df %>% 
  write_csv(here("pca_loadings_pk_coding.csv"))


##########

# calculate the variance for each gene
rv <- rowVars(assay(rld_filtered))

# Top n genes by variance to keep.
ntop <- 500

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- stats::prcomp(t(assay(rld_filtered)[select,]))

# Loadings for the first two PCs.
loadings <- pca$rotation[, seq_len(2)]

loadings_df <- as.data.frame(loadings) %>% 
  # mutate(gene_id = rownames(loadings))
  tibble::rownames_to_column(var = "gene_id") %>% 
  left_join(gene_annotation, by = "gene_id")

loadings_df %>% 
  write_csv(here("pca_loadings.csv"))
