library(here)

# read in (filtered) dataset
gene_counts_txi <- readRDS(here("results-refseq/r_objects/gene_counts_txi.RDS"))
gene_counts_txi_filtered <- readRDS(here("results-refseq/r_objects/gene_counts_txi_filtered.RDS"))

# extract Pk genes
gene_counts_txi_pk <- gene_counts_txi[gene_counts_txi$species == "pk",]
gene_counts_txi_filtered_pk <- gene_counts_txi_filtered[gene_counts_txi_filtered$species == "pk",]

# in case you need a clean table with only features and gene counts, you will need to drop the unnecessary columns
gene_counts_txi_clean <- subset(gene_counts_txi, select = -c(gene_name, gene_biotype, species, gene_type))
gene_counts_txi_pk_clean <- subset(gene_counts_txi_pk, select = -c(gene_name, gene_biotype, species, gene_type))
gene_counts_txi_filtered_clean <- subset(gene_counts_txi_filtered, select = -c(gene_name, gene_biotype, species, gene_type))
gene_counts_txi_filtered_pk_clean <- subset(gene_counts_txi_filtered_pk, select = -c(gene_name, gene_biotype, species, gene_type))

# to export the tables to csv
dir.create(here("results-refseq/cibersort"), )
write.csv(gene_counts_txi_clean, here("results-refseq/cibersort/gene_counts_txi_clean.csv"), row.names = FALSE)
write.csv(gene_counts_txi_pk_clean, here("results-refseq/cibersort/gene_counts_txi_pk_clean.csv"), row.names = FALSE)
write.csv(gene_counts_txi_filtered_clean, here("results-refseq/cibersort/gene_counts_txi_filtered_clean.csv"), row.names = FALSE)
write.csv(gene_counts_txi_filtered_pk_clean, here("results-refseq/cibersort/gene_counts_txi_filtered_pk_clean.csv"), row.names = FALSE)
