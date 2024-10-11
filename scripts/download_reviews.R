library(googlesheets4)
library(tidyverse)

gs4_auth("maxczap@ucsb.edu")

reviewers <- c(
  "tjl",
  "acn",
  "ahr",
  "arp",
  "cmh",
  "ecn",
  "jkb",
  "mfc",
  "tac"
)
links <- c(
"https://docs.google.com/spreadsheets/d/1r6GAoj-qsQvgRE6vQ-PRSvvbY7pUOVgIZtFVVBKG2Vs/edit?usp=drive_link",
"https://docs.google.com/spreadsheets/d/1JfLDnWzCt9W8yKWXQ3VmEDyHUutYmM4cVIb32bkLma0/edit?usp=drive_link",
"https://docs.google.com/spreadsheets/d/192zmhkjoRZ71d_2aDu2PappWWTM9f-CNeZ9j2TLIa5g/edit?usp=drive_link",
"https://docs.google.com/spreadsheets/d/1n5DuCEUHHWs8CUaaKrXAjBul6mUvol6JW5TF_xzTtHc/edit?usp=drive_link",
"https://docs.google.com/spreadsheets/d/1pLbDA5r9bDx2RoGoNg3pKJCKEGaOCFzkCNSDenSOY7g/edit?usp=drive_link",
"https://docs.google.com/spreadsheets/d/17TCcBig2jr4IrTAzRiRj6RifetFeLnlqEsDV5CW1Svk/edit?usp=drive_link",
"https://docs.google.com/spreadsheets/d/197tmebQiNloqYq2l_IXhs3YTnOsII9HwV6gE3QLRbC0/edit?usp=drive_link",
"https://docs.google.com/spreadsheets/d/1TGV72dx4_eGT-vJxEvBWvBElwpzAVdakfxMepDe3Lvs/edit?usp=drive_link",
"https://docs.google.com/spreadsheets/d/1PNfLK1fUxx17ufqDPppMQN-3B64NGu7um9zoKE9ny9w/edit?usp=drive_link"
)

# Reviews
reviews <- map2(links, reviewers,
                \(l, r) read_sheet(l,
                                   sheet = "papers",
                                   range = "a2:w8710",
                                   col_types = "c") %>%
                  mutate(source_sheet = r)) %>%
  list_rbind() %>%
  mutate(across(everything(), \(x) ifelse(x == "NA", NA, x))) %>%
  filter(assigned_to == source_sheet)
reviews$reviewed <- reviews[, 4:14] %>%
  lapply(is.na) %>%
  reduce(`&`) %>%
  `!`

reviews_file <- sprintf("reviews%s.rds", format(Sys.Date(), "%Y%m%d"))
saveRDS(reviews, file.path("outputs", "reviews", reviews_file))

# Additional taxa
add_taxa <- map(links,
                \(x) read_sheet(x,
                                sheet = "additional_taxa",
                                range = "a2:d1000",
                                col_types = "c")) %>%
  list_rbind() %>%
  drop_na()
add_taxa_file <- sprintf("taxa%s.rds", format(Sys.Date(), "%Y%m%d"))
saveRDS(add_taxa, file.path("outputs", "reviews", add_taxa_file))
