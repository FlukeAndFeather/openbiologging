library(tidyverse)

papers <- readxl::read_excel("data/systematicreview/scoring/reviewsample20230821mfc.xlsx",
                             sheet = "reviewsampletemplate20230821",
                             skip = 1) %>%
  filter(manuscript_type != "U",
         novel_biologging == "Y",
         biologging_context == "W") %>%
  select(id, authors:year)
devices <- readxl::read_excel("data/systematicreview/scoring/reviewsample20230821mfc.xlsx",
                              sheet = "devices",
                              skip = 1) %>%
  semi_join(papers, by = "id")
taxa <- readxl::read_excel("data/systematicreview/scoring/reviewsample20230821mfc.xlsx",
                           sheet = "taxa",
                           skip = 1) %>%
  semi_join(papers, by = "id")
data <- readxl::read_excel("data/systematicreview/scoring/reviewsample20230821mfc.xlsx",
                           sheet = "data",
                           skip = 1) %>%
  semi_join(papers, by = "id")

opendata <- papers %>%
  right_join(data, by = "id") %>%
  left_join(taxa, by = "id", relationship = "many-to-many") %>%
  left_join(devices, by = "id")

fair_df <- opendata %>%
  # If no availability statement, then not FAIR
  pivot_longer(F1:"R1.3", names_to = "fair_attr", values_to = "value") %>%
  mutate(value = ifelse(availability_statement == "N", "N", value)) %>%
  pivot_wider(names_from = fair_attr, values_from = value) %>%
  mutate(year_2007 = year - 2007,
         across(c(is_biologging, F1:"R1.3"),
                \(x) case_match(x,
                                "Y" ~ TRUE,
                                "N" ~ FALSE,
                                .default = NA))) %>%
  drop_na(F1:"R1.3")

fair_long <- fair_df %>%
  pivot_longer(F1:"R1.3", names_to = "fair_attr", values_to = "value") %>%
  mutate(is_biologging = ifelse(is_biologging, "Biologging", "Other"))

ggplot(fair_long,
       aes(x = interaction(id, description), y = fair_attr, color = value)) +
  geom_point() +
  scale_color_manual("Meets requirement",
                     values = c("TRUE" = "cornflowerblue",
                                "FALSE" = "firebrick")) +
  facet_grid(rows = "is_biologging") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

# I think we need to use GEEs to account for correlations between FAIR attributes

fair_scores <- fair_long %>%
  group_by(id, description) %>%
  summarize(fair_attrs = sum(value), .groups = "drop",
            fair_n = n()) %>%
  right_join(distinct(select(fair_df, id, description, is_biologging, year)),
             by = c("id", "description")) %>%
  mutate(is_biologging = ifelse(is_biologging, "Biologging", "Other"))
ggplot(fair_scores, aes(year, fair_attrs / fair_n)) +
  geom_point(position = "jitter") +
  geom_smooth(method = "gam", formula = y ~ x) +
  facet_grid(rows = "is_biologging") +
  theme_minimal()
fair_model <- lme4::glmer(value ~ year_2007 * is_biologging + (1 | fair_attr),
                          fair_long,
                          family = "binomial")
summary(fair_model)
emmeans::emtrends(fair_model, ~ is_biologging, "year_2007") %>%
  plot()
