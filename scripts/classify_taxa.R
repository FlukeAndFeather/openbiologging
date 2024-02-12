library(taxize)
library(tidyverse)

# Read reviews and additional taxa
latest_reviews <- dir("outputs/reviews", "reviews.*rds", full.names = TRUE) %>%
  sort(decreasing = TRUE) %>%
  first() %>%
  readRDS()
latest_add_taxa <- dir("outputs/reviews", "taxa.*rds", full.names = TRUE) %>%
  sort(decreasing = TRUE) %>%
  first() %>%
  readRDS()

novel_biologging <- latest_reviews %>%
  filter(manuscript_type != "U",
         novel_biologging == "Y",
         str_detect(biologging_context, "W"),
         reviewed)

biologging_taxa <- select(novel_biologging, genus, species) %>%
  rbind(latest_add_taxa %>%
          filter(id %in% novel_biologging$id) %>%
          select(genus, species)) %>%
  distinct() %>%
  arrange(genus, species) %>%
  mutate(sciname = paste(genus, species))

# Write taxa to text file

writeLines(c("name", biologging_taxa$sciname),
           "outputs/taxa/biologgingtaxa.txt")

# Use taxmatch to get tsn's https://itis.gov/taxmatch.html

# Read matches
taxa_matches <- read_delim("outputs/taxa/matchedtaxa.txt", delim = "|") %>%
  drop_na(TSN) %>%
  mutate(accepted_tsn = coalesce(`Accepted TSN`, TSN))

# Get classifications
taxa_class <- classification(taxa_matches$accepted_tsn, db = "itis")
taxa_class_df <- map(
  taxa_class,
  \(x) {
    x %>%
      filter(rank %in% c("phylum", "class", "order", "family", "species")) %>%
      select(rank, name) %>%
      pivot_wider(names_from = "rank", values_from = "name")
  }
) %>%
  list_rbind() %>%
  mutate(id = names(taxa_class)) %>%
  relocate(id) %>%
  distinct() %>%
  right_join(transmute(taxa_matches,
                       id = as.character(accepted_tsn),
                       orig_id = TSN,
                       query = `Scientific Name`),
             by = "id",
             relationship = "one-to-many") %>%
  select(-orig_id) %>%
  distinct()

# Join classifications with studies
study_taxa <- rbind(
  select(novel_biologging, id, genus, species),
  select(filter(latest_add_taxa, id %in% novel_biologging$id), id, genus, species)
) %>%
  transmute(id, sciname = paste(genus, species)) %>%
  left_join(select(taxa_class_df, -id),
            by = c(sciname = "query"),
            relationship = "many-to-one") %>%
  drop_na(phylum:species)
saveRDS(study_taxa, "outputs/taxa/study_taxa.rds")
