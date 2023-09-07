library(googlesheets4)
library(tidyverse)


# Read data ---------------------------------------------------------------
gs4_deauth() # Just so Google doesn't try to authenticate you}
read_review <- function(sheet) {
  read_sheet(
    "https://docs.google.com/spreadsheets/d/1zRA77pWp4nFZ9sjhVKXxxe0TBSg96QBC0L8XThshfQE/edit?usp=sharing",
    sheet = sheet,
    skip = 1,
    na = c("", "NA")
  )
}
papers <- read_review("papers")
devices <- read_review("devices")
taxa <- read_review("taxa")
datasources <- read_review("data")


# Clean data --------------------------------------------------------------
novel_data_papers <- papers %>%
  drop_na(manuscript_type) %>%
  filter(manuscript_type != "U",
         novel_biologging == "Y",
         biologging_context == "W")
novel_data_devices <- devices %>%
  select(-note) %>%
  semi_join(novel_data_papers, by = "id")
novel_data_taxa <- taxa %>%
  select(-note) %>%
  semi_join(novel_data_papers, by = "id")
novel_data_sources <- datasources %>%
  select(-note) %>%
  semi_join(novel_data_papers, by = "id") %>%
  mutate(across(c(is_biologging:availability_statement, F1:Complete),
                \(x) case_match(x, "Y" ~ TRUE, "N" ~ FALSE, .default = NA)))

# Save outputs ------------------------------------------------------------
saveRDS(papers, "outputs/papers.rds")
saveRDS(devices, "outputs/devices.rds")
saveRDS(taxa, "outputs/taxa.rds")
saveRDS(datasources, "outputs/datasources.rds")
saveRDS(novel_data_papers, "outputs/novel_data_papers.rds")
saveRDS(novel_data_devices, "outputs/novel_data_devices.rds")
saveRDS(novel_data_taxa, "outputs/novel_data_taxa.rds")
saveRDS(novel_data_sources, "outputs/novel_data_sources.rds")
