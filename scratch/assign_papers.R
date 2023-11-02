library(tidyverse)
library(googlesheets4)

# gs4_auth("mczapans@ucsc.edu")

all_papers <- read_sheet("https://docs.google.com/spreadsheets/d/1C1P-wCDw29eZBYwn-i1njMssvQrF6s0unZuJaXn8gls")

collaborators <- c("arp", "ecn", "cmh", "jkb", "mfc", "star", "tac")

set.seed(2007)
assignments <- all_papers %>%
  rbind(all_papers) %>%
  mutate(assigned_to = sample(collaborators, size = nrow(.), replace = TRUE)) %>%
  arrange(assigned_to)

write_sheet(assignments, "https://docs.google.com/spreadsheets/d/1432iil-JwVyY2l1BqH3tKAlRKKKZlKdBxueSQy-oCO4")
