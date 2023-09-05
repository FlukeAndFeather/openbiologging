library(taxize)

taxa <- readRDS("outputs/taxa.rds")
# Make sure ENTREZ_KEY is set in .Renviron; see taxize::use_entrez()
sci_names <- paste(taxa$genus, taxa$species) %>%
  unique() %>%
  sort()
sci_names <- sci_names[!sci_names %in% c("??? ???", "NA NA")]
t <- tax_name(sci = sci_names,
              get = c("phylum", "class", "order", "family", "genus"),
              db = "itis") %>%
  mutate(species = str_extract(query, "[A-Z][a-z]+ ([a-z]+)", group = 1))
saveRDS(t, "outputs/taxa_classification.rds")
