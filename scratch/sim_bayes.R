library(brms)
library(tidyverse)

set.seed(1234)

# Grid of independent variables
biolog_grid <- expand_grid(
  era = c("early", "middle", "recent"),
  E = c("marine", "terrestrial"),
  S = c("aspatial", "spatial")
) %>%
  slice(rep(1:12, each = 50)) %>%
  mutate(T = case_when(
    era == "early"  ~ sample(2007:2015, 600, replace = TRUE),
    era == "middle" ~ sample(2016:2019, 600, replace = TRUE),
    era == "recent" ~ sample(2020:2023, 600, replace = TRUE)
  ))

# Parameters for simulation
# delta
delta <- c(
  # Dummy for 2007
  0,
  # slow increases in early period
  rep(1, 8),
  # fastest increases in middle period
  rep(3, 4),
  # moderate increases in recent period
  rep(2, 4)
)
delta <- delta / sum(delta)
# betas (E, S, T)
beta_E <- c(marine = 0, terrestrial = 1)
beta_S <- c(aspatial = 0, spatial = 1)
beta_T <- 3
alpha <- -4

# Run simulation
inv_logit <- \(x) exp(x) / (1 + exp(x))
biolog <- biolog_grid %>%
  mutate(cumsum_delta = map_dbl(T, \(t) sum(delta[1:(t - 2007)])),
         phi = alpha + beta_E[E] + beta_S[S] + beta_T * cumsum_delta,
         F = rbinom(600, size = 1, prob = inv_logit(phi)))

# Visualize data
biolog_summ <- biolog %>%
  group_by(T, E, S) %>%
  summarize(F = mean(F),
            .groups = "drop")
# phi ~ T
ggplot(biolog, aes(T, phi, color = E)) +
  geom_point() +
  facet_wrap(~S) +
  theme_classic() +
  theme(legend.position = "bottom")
# F ~ T
ggplot(biolog, aes(T, F, color = E)) +
  geom_jitter(aes(shape = S), alpha = 0.75, width = 0.25, height = 0.05) +
  geom_line(aes(linetype = S), biolog_summ) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_classic() +
  theme(legend.position = "bottom")

# Fit model
biolog_prior <- c(
  set_prior(prior = "normal(0, 1)", class = "Intercept"),
  set_prior(prior = "normal(0, 1)", coef = "Eterrestrial"),
  set_prior(prior = "normal(0, 1)", coef = "Sspatial"),
  set_prior(prior = "normal(0, 1)", coef = "moT"),
  set_prior(prior = "dirichlet(1)", class = "simo", coef = "moT1")
)
biolog_mod <- brm(F ~ E + S + mo(T),
                  data = biolog,
                  family = bernoulli(link = "logit"),
                  prior = biolog_prior)
posterior_summary(biolog_mod) %>%
  round(digits = 3)

# Hypothesis testing
## H1
### Calculate contrasts
h1_data_2016 <- expand_grid(
  E = c("marine", "terrestrial"),
  S = c("aspatial", "spatial"),
  T = 2016
)
h1_pred_2016 <- posterior_predict(biolog_mod,
                                  newdata = h1_data_2016,
                                  ndraws = 1000)
h1_data_2023 <- expand_grid(
  E = c("marine", "terrestrial"),
  S = c("aspatial", "spatial"),
  T = 2023
)
h1_pred_2023 <- posterior_predict(biolog_mod,
                                  newdata = h1_data_2023,
                                  ndraws = 1000)
h1_contrast <- (h1_pred_2023 - h1_pred_2016) %>%
  t() %>%
  cbind(select(h1_data_2016, E, S), .) %>%
  pivot_longer(-c(E, S), names_to = "draw", values_to = "contrast")
### Plot contrasts
ggplot(h1_contrast, aes(contrast, fill = S)) +
  geom_bar(aes(y = after_stat(prop)), position = "dodge") +
  facet_wrap(~E) +
  scale_fill_manual(values = c("firebrick", "cornflowerblue")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Contrast of F\n(2023 vs 2016)",
       y = "Proportion",
       fill = "Sensor type") +
  theme_classic()
### Table of contrasts
h1_contrast %>%
  count(E, S, contrast) %>%
  group_by(E, S) %>%
  mutate(n = n / sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = "contrast", values_from = "n") %>%
  mutate(across(-c(E, S), scales::percent))

## H2
### Calculate contrasts
h2_data_marine <- expand_grid(
  E = "marine",
  S = c("aspatial", "spatial"),
  T = 2023
)
h2_pred_marine <- posterior_predict(biolog_mod,
                                    newdata = h2_data_marine,
                                    ndraws = 1000)
h2_data_terrestrial <- expand_grid(
  E = "terrestrial",
  S = c("aspatial", "spatial"),
  T = 2023
)
h2_pred_terrestrial <- posterior_predict(biolog_mod,
                                         newdata = h2_data_terrestrial,
                                         ndraws = 1000)
h2_contrast <- (h2_pred_terrestrial - h2_pred_marine) %>%
  t() %>%
  cbind(select(h2_data_marine, S), .) %>%
  pivot_longer(-c(S), names_to = "draw", values_to = "contrast")
### Plot contrasts
ggplot(h2_contrast, aes(contrast, fill = S)) +
  geom_bar(aes(y = after_stat(prop)), position = "dodge") +
  scale_fill_manual(values = c("firebrick", "cornflowerblue")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Contrast of F\n(Terrestrial vs marine)",
       y = "Proportion",
       fill = "Sensor type") +
  theme_classic()
### Table of contrasts
h2_contrast %>%
  count(S, contrast) %>%
  group_by(S) %>%
  mutate(n = n / sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = "contrast", values_from = "n") %>%
  mutate(across(-S, scales::percent))

## H3
### Calculate contrasts
h3_data_aspatial <- expand_grid(
  E = c("marine", "terrestrial"),
  S = "aspatial",
  T = 2023
)
h3_pred_aspatial <- posterior_predict(biolog_mod,
                                      newdata = h3_data_aspatial,
                                      ndraws = 1000)
h3_data_spatial <- expand_grid(
  E = c("marine", "terrestrial"),
  S = "spatial",
  T = 2023
)
h3_pred_spatial <- posterior_predict(biolog_mod,
                                     newdata = h3_data_spatial,
                                     ndraws = 1000)
h3_contrast <- (h3_pred_spatial - h3_pred_aspatial) %>%
  t() %>%
  cbind(select(h3_data_aspatial, E), .) %>%
  pivot_longer(-E, names_to = "draw", values_to = "contrast")
### Plot contrasts
ggplot(h3_contrast, aes(contrast, fill = E)) +
  geom_bar(aes(y = after_stat(prop)), position = "dodge") +
  scale_fill_manual(values = c("mediumpurple", "tan")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Contrast of F\n(Spatial vs aspatial)",
       y = "Proportion",
       fill = "Ecosystem") +
  theme_classic()
### Table of contrasts
h3_contrast %>%
  count(E, contrast) %>%
  group_by(E) %>%
  mutate(n = n / sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = "contrast", values_from = "n") %>%
  mutate(across(-E, scales::percent))
