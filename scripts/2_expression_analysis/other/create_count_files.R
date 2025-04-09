library(here)
library(dplyr)
library(readr)
library(stringr)
library(tximport)
library(DESeq2)

gene_tpm_filtered_long %>%
  left_join(samples, by = "sample") %>%
  group_by(sample_name, gene_type) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads)) %>% 
  write_csv(here("./tpm_counts_general.csv"))

gene_tpm_filtered_long %>%
  mutate(gene_biotype = coalesce(gene_biotype, gene_type)) %>% 
  left_join(samples, by = "sample") %>%
  group_by(sample_name, gene_biotype) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads)) %>% 
  write_csv(here("./tpm_counts_genebiotypes.csv"))

gene_tpm_filtered_long %>%
  left_join(samples, by = "sample") %>%
  filter(species == "pk" & gene_biotype == "protein_coding") %>%
  group_by(sample_name) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads)) %>% 
  write_csv(here("./tpm_counts_pkcoding.csv"))



gene_counts_txi_filtered_long %>%
  left_join(samples, by = "sample") %>%
  group_by(sample_name, gene_type) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads)) %>% 
  write_csv(here("./raw_counts_general.csv"))

gene_counts_txi_filtered_long %>%
  mutate(gene_biotype = coalesce(gene_biotype, gene_type)) %>% 
  left_join(samples, by = "sample") %>%
  group_by(sample_name, gene_biotype) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads)) %>% 
  write_csv(here("./raw_counts_genebiotypes.csv"))

gene_counts_txi_filtered_long %>%
  left_join(samples, by = "sample") %>%
  filter(species == "pk" & gene_biotype == "protein_coding") %>%
  group_by(sample_name) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads)) %>% 
  write_csv(here("./raw_counts_pkcoding.csv"))

###

# data_list <- list.files(path = here("results/r_objects"),
#                         pattern = ".RDS",
#                         recursive = TRUE, full.names = TRUE) %>%
#   purrr::set_names(., stringr::str_remove(basename(.), ".RDS")) %>%
#   purrr::imap(~ readRDS(file = .x ))

# salmon output files

gene_counts_txi <- read_delim(here("results/star_salmon/salmon.merged.gene_counts.tsv")) %>%
  mutate(across(starts_with("X"), \(x) as.integer(round(x)))) %>%
  rename_with(~ str_replace_all(.x, pattern = "\\.", replacement = "-"), .cols = `X104974.001.001_S1_L001`:`X104974.001.016_S1_L001`) %>%
  rename_with(~ str_sub(.x, 2, -1), .cols = `X104974-001-001_S1_L001`:`X104974-001-016_S1_L001`)

transcript_counts_txi <- read_delim(here("results/star_salmon/salmon.merged.transcript_counts.tsv")) %>%
  mutate(across(starts_with("X"), \(x) as.integer(round(x)))) %>%
  rename_with(~ str_replace_all(.x, pattern = "\\.", replacement = "-"), .cols = `X104974.001.001_S1_L001`:`X104974.001.016_S1_L001`) %>%
  rename_with(~ str_sub(.x, 2, -1), .cols = `X104974-001-001_S1_L001`:`X104974-001-016_S1_L001`)

gene_counts_txi %>% summarise(across(starts_with("104974"), ~ sum(., na.rm = TRUE)))
transcript_counts_txi %>% summarise(across(starts_with("104974"), ~ sum(., na.rm = TRUE)))


# tximport approach

salmon_quant_files <- paste0("./results/star_salmon/", list.files(path = "./results/star_salmon/", pattern = "quant.sf", recursive = TRUE))
transcript_info <- read.csv(here("results/star_salmon/tx2gene.tsv"),
                            sep = "\t", header = FALSE,
                            col.names = c("tx", "gene_id", "gene_name")
)
transcript_info <- list(
  transcript = transcript_info,
  gene = unique(transcript_info[, 2:3]),
  tx2gene = transcript_info[, 1:2]
)
txi <- tximport(salmon_quant_files, type = "salmon", tx2gene = transcript_info$tx2gene, txOut = FALSE)
colSums(txi$counts)

samples <- read_csv(here("./data/samplesheet.csv"),
                    col_select = c("sample", "condition", "RNA_type", "library_kit", "removal", "WBC_depletion", "kit"),
                    col_types = cols(.default = "f")
) %>%
  mutate(RNA_type = str_to_title(RNA_type)) %>%
  mutate(filter_name = factor(WBC_depletion,
                              levels = c("nofilter", "centpip", "plasmo", "pmacs", "cellulose"),
                              labels = c("no filter", "cent/pip", "plasmodipur", "PMACS", "cellulose")
  )) %>%
  # mutate(kit = as.factor(str_replace_all(kit, "_", " "))) %>%
  mutate(kit_name = factor(kit,
                           levels = c("mRNA_qiagen_FSglob", "mRNA_qiagen_FSglobH", "tRNA_qiagen_FSglobH", "mRNA_illumina"),
                           labels = c("mRNA Qiagen FS globin", "mRNA Qiagen FS globin/rRNA", "tRNA Qiagen FS globin/rRNA", "mRNA Illumina")
  )) %>%
  mutate(sample_name = as.factor(str_replace(str_replace(sample, "104974-001-0", "Sample "), "_S1_L001", paste(" -", kit_name, "-", filter_name)))) %>%
  mutate(sample_short = as.factor(str_replace(str_replace(sample, "104974-001-0", "Sample "), "_S1_L001", "")))
gse <- DESeqDataSetFromTximport(txi, colData = samples, design = ~1)
colSums(assay(gse))




data_list <- list.files(path = here("results/r_objects"),
                        pattern = ".RDS",
                        recursive = TRUE, full.names = TRUE) %>%
  purrr::set_names(., stringr::str_remove(basename(.), ".RDS")) %>%
  purrr::imap(~ readRDS(file = .x ))

data_list$gene_counts_txi_filtered_long %>%
  left_join(samples, by = "sample") %>%
  # group_by(sample_name, gene_type) %>%
  group_by(sample_name) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads))

data_list$gene_counts_txi_long %>%
  left_join(samples, by = "sample") %>%
  group_by(sample_name) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads))




data_list$gene_counts_txi_filtered_long %>%
  left_join(samples, by = "sample") %>%
  filter(species == "pk") %>%
  # group_by(sample_name, gene_type) %>%
  group_by(sample_name) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads))

data_list$gene_counts_txi_long %>%
  left_join(samples, by = "sample") %>%
  filter(species == "pk") %>%
  # group_by(sample_name, gene_type) %>%
  group_by(sample_name) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads))


data_list$gene_counts_txi_filtered_long %>%
  left_join(samples, by = "sample") %>%
  filter(species == "pk" & gene_biotype == "protein_coding") %>%
  # group_by(sample_name, gene_type) %>%
  group_by(sample_name) %>%
  summarise(sum_reads = sum(counts, na.rm = TRUE)) %>%
  mutate(per = 100 * sum_reads / sum(sum_reads))


# transcript counts for pk genes

pk_genes <- data_list$gene_counts_txi_filtered %>% 
  filter(species == "pk") %>% 
  select(gene_id,gene_name)

transcript_counts_txi %>% inner_join(pk_genes, by="gene_id") %>% 
  summarise(across(starts_with("104974"), ~ sum(., na.rm = TRUE)))

# gene counts for pk

gene_counts_txi %>% inner_join(pk_genes, by="gene_id") %>% 
  summarise(across(starts_with("104974"), ~ sum(., na.rm = TRUE)))
colSums(gene_counts_txi %>% inner_join(pk_genes, by="gene_id") %>% select(starts_with("104974")))
