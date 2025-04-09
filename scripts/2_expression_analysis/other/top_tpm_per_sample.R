a<-gene_tpm_filtered_long %>% 
  filter(species == "pk") %>% 
  slice_max(n = 20, order_by = counts, by = sample) %>% 
  write_csv(here("top_10_tpm_per_sample.csv"))

gene_tpm_filtered_long %>% 
  filter(species == "pk") %>% 
  group_by(sample) %>% 
  top_n(wt = counts, n = 10) %>% 
  select(c(gene_name, sample, counts))
