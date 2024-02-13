library(brms)
library(tidyverse)

set.seed(702)

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

# Contrasts
pred_grid <- expand_grid(
  E = c("marine", "terrestrial"),
  S = c("aspatial", "spatial"),
  T = 2007:2023
)
biolog_pred_wide <- posterior_predict(biolog_mod,
                                      newdata = pred_grid,
                                      ndraws = 1000)
biolog_pred <- cbind(pred_grid, t(biolog_pred_wide)) %>%
  pivot_longer(-c(E, S, T), names_to = "draw", values_to = "F")
ggplot(biolog_pred, aes(T, F, color = E)) +
  geom_line(aes(linetype = S), biolog_summ) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  theme_classic() +
  theme(legend.position = "bottom")

# Ecosystem: marine and terrestrial (in 2023)
E_contrast <- biolog_pred %>%
  filter(T == 2023) %>%
  group_by(draw, S) %>%
  summarize(contrast = F[E == "terrestrial"] - F[E == "marine"],
            .groups = "drop")

ggplot(E_contrast, aes(contrast, fill = S)) +
  geom_bar(position = "dodge") +
  labs(x = "Contrast of F (presence of permanent identifier)\n(Terrestrial relative to marine)",
       y = "Count\n(1000 draws from posterior distribution)",
       fill = "Sensor type") +
  theme_classic()

E_contrast %>%
  group_by(S) %>%
  summarize(`Mean contrast` = mean(contrast))

# Sensor: aspatial and spatial
S_contrast <- biolog_pred %>%
  filter(T == 2023) %>%
  group_by(draw, E) %>%
  summarize(contrast = F[S == "spatial"] - F[S == "aspatial"],
            .groups = "drop")

ggplot(S_contrast, aes(contrast, fill = E)) +
  geom_bar(position = "dodge") +
  labs(x = "Contrast of F (presence of permanent identifier)\n(Spatial relative to aspatial)",
       y = "Count\n(1000 draws from posterior distribution)",
       fill = "Ecosystem") +
  theme_classic()

S_contrast %>%
  group_by(E) %>%
  summarize(`Mean contrast` = mean(contrast))

# Time: 2023 vs 2013
T_contrast <- biolog_pred %>%
  filter(T %in% c(2013, 2023)) %>%
  group_by(draw, E, S) %>%
  summarize(contrast = mean(F[T == 2023] - F[T == 2013]),
            .groups = "drop")

ggplot(T_contrast, aes(contrast, fill = E)) +
  geom_bar(position = "dodge") +
  facet_wrap(~S) +
  labs(x = "Contrast of F (presence of permanent identifier)\n(2023 relative to 2013)",
       y = "Count\n(1000 draws from posterior distribution)",
       fill = "Ecosystem") +
  theme_classic()

T_contrast %>%
  group_by(E, S) %>%
  summarize(`Mean contrast` = mean(contrast),
            .groups = "drop") %>%
  pivot_wider(names_from = "S", values_from = "Mean contrast")
