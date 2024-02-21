library(tidyverse)

# hypotheses --------------------------------------------------------------

tibble(Time = 1:20,
       Open = cumsum(rep(c(1, -1, 1, 1), 5))) %>%
  mutate(Open = Open / max(Open) * 0.8) %>%
  ggplot(aes(Time, Open)) +
  geom_line(arrow = arrow(), linewidth = 2) +
  scale_y_continuous("Studies with open data",
                     labels = scales::percent,
                     limits = c(0, 1)) +
  theme_classic(18) +
  theme(axis.text.x = element_blank())
ggsave("~/Downloads/bls8_h1.svg", width = 126.5, height = 90, units = "mm")

tibble(Ecosystem = c("Marine", "Terrestrial"),
       Open = c(0.4, 0.7),
       Err = c(0.08, 0.15)) %>%
  ggplot(aes(Ecosystem, Open, fill = Ecosystem)) +
  geom_col() +
  geom_pointrange(aes(ymin = Open - Err, ymax = Open + Err),
                  size = 1.5,
                  linewidth = 2) +
  scale_y_continuous("Studies with open data",
                     labels = scales::percent,
                     limits = c(0, 1)) +
  scale_fill_manual(values = c("#3070AD", "#469C76")) +
  theme_classic(18) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")
ggsave("~/Downloads/bls8_h2.svg", width = 126.5, height = 90, units = "mm")

tibble(Sensor = c("Aspatial", "Spatial"),
       Open = c(0.4, 0.7),
       Err = c(0.08, 0.15)) %>%
  ggplot(aes(Sensor, Open, fill = Sensor)) +
  geom_col() +
  geom_pointrange(aes(ymin = Open - Err, ymax = Open + Err),
                  size = 1.5,
                  linewidth = 2) +
  scale_y_continuous("Studies with open data",
                     labels = scales::percent,
                     limits = c(0, 1)) +
  scale_fill_manual(values = c("#DE6E79", "#C9BC59")) +
  theme_classic(18) +
  theme(axis.title.x = element_blank(),
        legend.position = "none")
ggsave("~/Downloads/bls8_h3.svg", width = 126.5, height = 90, units = "mm")

# results -----------------------------------------------------------------

library(tidyverse)

# Get the most recent set of reviews
reviews <- dir(here::here("outputs", "reviews"),
               "reviews.*rds",
               full.names = TRUE) %>%
  sort(decreasing = TRUE) %>%
  first() %>%
  readRDS()

# And most recent studies-by-taxonomy
study_taxa <- readRDS(here::here("outputs", "taxa", "study_taxa.rds"))

scored_papers <- reviews %>%
  mutate(across(everything(), \(x) ifelse(x == "NA", NA, x))) %>%
  filter(reviewed)

nrow(scored_papers)

scored_summary <- scored_papers %>%
  mutate(related = manuscript_type != "U" &
           str_detect(biologging_context, "W"),
         novel = related & novel_biologging == "Y",
         das = related & biologging_availability == "Y") %>%
  summarize(
    total = n(),
    related = sum(related, na.rm = TRUE),
    novel = sum(novel, na.rm = TRUE),
    das = sum(das, na.rm = TRUE))

scored_papers %>%
  summarize(
    `Total papers` = n(),
    `Related` = sum(manuscript_type != "U", na.rm = TRUE),
    `Novel biologging` = sum(novel_biologging == "Y",
                             na.rm = TRUE),
    `Data availability` = sum(biologging_availability == "Y",
                              na.rm = TRUE)) %>%
  pivot_longer(everything()) %>%
  mutate(name = fct_reorder(name, value, .fun = \(x) -x)) %>%
  ggplot(aes(name, value)) +
  geom_col() +
  labs(y = "Papers") +
  theme_classic() +
  theme(axis.title.x = element_blank())

novel_biolog <- scored_papers %>%
  filter(manuscript_type != "U",
         novel_biologging == "Y",
         str_detect(biologging_context, "W"),
         !is.na(novel_biologging))

novel_biolog %>%
  mutate(Habitat = case_when(
    substr(habitat, 1, 1) == "M" ~ "Marine",
    substr(habitat, 1, 1) == "T" ~ "Terrestrial",
    TRUE ~ "Other"
  )) %>%
  group_by(Habitat) %>%
  summarize(
    Overall = n(),
    Location = sum(str_detect(device_cat, "L"),
                   na.rm = TRUE),
    Intrinsic = sum(str_detect(device_cat, "I"),
                    na.rm = TRUE),
    Environment = sum(str_detect(device_cat, "E"),
                      na.rm = TRUE)) %>%
  pivot_longer(-Habitat,
               names_to = "Sensor type",
               values_to = "Papers") %>%
  mutate(
    `Sensor type` = fct_reorder(`Sensor type`,
                                Papers,
                                .fun = \(x) -sum(x)),
    Habitat = fct_reorder(Habitat, Papers, .fun = sum)
  ) %>%
  ggplot(aes(`Sensor type`, Papers, fill = Habitat)) +
  geom_col() +
  scale_fill_manual(values = c(Terrestrial = "tan1",
                               Marine = "cornflowerblue",
                               Other = "darkorchid2")) +
  theme_classic(base_size = 24) +
  theme(legend.position = c(0.9, 0.9),
        legend.justification = c(1, 1))
ggsave("scratch/sensors.svg", width = 253, height = 213, units = "mm")

data_avail <- novel_biolog %>%
  drop_na(biologging_availability) %>%
  mutate(
    biologging_availability = ifelse(biologging_availability == "Y",
                                     1, 0),
    year = as.numeric(year),
    year2007 = year - 2007
  )

data_avail_mod <- glm(biologging_availability ~ year2007,
                      family = "binomial",
                      data = data_avail)
data_avail_grid <- tibble(year2007 = 0:16)
data_avail_pred <- predict(data_avail_mod,
                           newdata = data_avail_grid,
                           se.fit = TRUE)
invlogit <- binomial()$linkinv
data_avail_pred_df <- data_avail_grid %>%
  mutate(eta = data_avail_pred$fit,
         eta_lwr = data_avail_pred$fit - data_avail_pred$se.fit,
         eta_upr = data_avail_pred$fit + data_avail_pred$se.fit,
         biologging_availability = invlogit(eta),
         biologging_availability_lwr = invlogit(eta_lwr),
         biologging_availability_upr = invlogit(eta_upr),
         year = year2007 + 2007)

ggplot(data_avail, aes(year, biologging_availability)) +
  geom_point(shape = 21,
             position = position_jitter(width = 0.2, height = 0.05)) +
  geom_ribbon(aes(x = year,
                  ymin = biologging_availability_lwr,
                  ymax = biologging_availability_upr),
              data_avail_pred_df,
              alpha = 0.5) +
  geom_line(data = data_avail_pred_df,
            color = "blue") +
  labs(x = "Year", y = "Data Availability Statement") +
  theme_classic(base_size = 24)
ggsave("scratch/das.svg", width = 253, height = 213, units = "mm")


study_taxa_categorized <- study_taxa %>%
  mutate(group1 = ifelse(phylum == "Chordata", "Vertebrate", "Invertebrate"),
         group2 = case_when(
           class == "Aves" ~ "Birds",
           class == "Mammalia" ~ "Mammals",
           class %in% c("Chondrichthyes", "Chondrostei", "Teleostei") ~ "Fish",
           group1 == "Vertebrate" ~ "Other vertebrates",
           TRUE ~ phylum
         ),
         group3 = case_when(
           order == "Carnivora" ~ "Carnivores",
           order %in% c("Artiodactyla", "Perissodactyla") ~ "Ungulates",
           order == "Cetacea" ~ "Cetaceans",
           class == "Mammalia" ~ "Other mammals",
           order %in% c("Charadriiformes", "Sphenisciformes", "Procellariiformes", "Suliformes", "Pelecaniformes") ~ "Seabirds",
           order %in% c("Anseriformes", "Galliformes") ~ "Fowl",
           order %in% c("Accipitriformes", "Falconiformes", "Strigiformes") ~ "Raptors",
           class == "Aves" ~ "Other birds",
           TRUE ~ group2
         ),
         group2_fill = group2) %>%
  to_lodes_form(group1:group3)

ggplot(study_taxa_categorized, aes(x = x, stratum = stratum, alluvium = alluvium)) +
  stat_alluvium(aes(fill = group2_fill), decreasing = FALSE) +
  stat_stratum(decreasing = FALSE) +
  stat_stratum(aes(label = sprintf("%s (%d)", stratum, after_stat(n))),
                   geom = "text", decreasing = FALSE, size = 6) +
  theme_void() +
  theme(legend.position = "none")
ggsave("scratch/sankey.svg", width = 253, height = 202, units = "mm")



