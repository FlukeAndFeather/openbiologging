library(tidyverse)
library(googlesheets4)
gs4_auth("mczapans@ucsc.edu")

all_papers <- read_sheet("https://docs.google.com/spreadsheets/d/1C1P-wCDw29eZBYwn-i1njMssvQrF6s0unZuJaXn8gls") %>%
  mutate(volume = as.double(volume),
         issue = as.double(issue),
         id = as.character(id))
already_reviewed <- read_csv("data/systematicreview/scoring/ReviewsInProgress.csv")
# some papers already reviewed by 3 individuals; drop mfc and keep others
already_reviewed %>%
  group_by(paperid) %>%
  summarize(n = n(),
            reviewers = paste(reviewer_initials, collapse = ", ")) %>%
  arrange(desc(n))
already_reviewed <- already_reviewed %>%
  anti_join(already_reviewed %>%
              count(paperid) %>%
              filter(n == 3) %>%
              mutate(reviewer_initials = "mfc"),
            by = c("paperid", "reviewer_initials")) %>%
  group_by(paperid) %>%
  mutate(nth = seq(n())) %>%
  ungroup()
already_reviewed %>%
  group_by(paperid) %>%
  summarize(n = n(),
            reviewers = paste(reviewer_initials, collapse = ", ")) %>%
  arrange(desc(n))

collaborators <- c("acn", "arp", "ecn", "cmh", "jkb", "mfc", "star", "tac")

set.seed(2007)
assignments <- mutate(all_papers, nth = 1) %>%
  rbind(mutate(all_papers, nth = 2)) %>%
  left_join(select(already_reviewed, id = paperid, reviewer_initials, nth),
            by = c("id", "nth")) %>%
  mutate(
    assigned_to = coalesce(
      reviewer_initials,
      sample(collaborators, size = nrow(.), replace = TRUE)
    )
  ) %>%
  arrange(assigned_to)

write_sheet(assignments, "https://docs.google.com/spreadsheets/d/1432iil-JwVyY2l1BqH3tKAlRKKKZlKdBxueSQy-oCO4")
