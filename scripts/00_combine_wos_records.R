library(readxl)
library(tidyverse)

wos_dir <- "data/systematicreview/wosexport20230821/"
make_id <- function(authors, year) {
  aut <- author %>%
    str_extract("[^,]*") %>%
    str_to_lower() %>%
    str_replace_all("[^a-zA-Z]", "")
  tibble(id = paste0(aut, year)) %>%
    group_by(id) %>%
    mutate(id = if(n() > 1) paste0(id, letters[1:n()]) else id) %>%
    pull(id)
}
wos_records <- dir(wos_dir, full.names = TRUE) %>%
  map(read_excel, col_types = "text") %>%
  list_rbind() %>%
  transmute(manuscript_type = "",
            novel_biologging = "",
            other_data = "",
            biologging_fair = "",
            other_fair = "",
            authors = Authors,
            title = `Article Title`,
            source = `Source Title`,
            year = `Publication Year`,
            volume = Volume,
            issue = Issue,
            doi = DOI,
            pub_type = `Publication Type`,
            id = make_id(authors, year))
write_csv(wos_records, file.path(wos_dir, "reviewtemplate20230821.csv"))
