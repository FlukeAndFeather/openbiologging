library(brms)
library(cowplot)
library(tidybayes)
library(tidyverse)
set.seed(1234)

# Grid of independent variables
biolog_grid <- expand_grid(
  era = c("early", "middle", "recent"),
  E = c("marine", "terrestrial"),
  S = c("aspatial", "spatial")
) %>%
  slice(rep(seq(nrow(.)), each = 50)) %>%
  mutate(.id = seq(nrow(.)),
         T = case_when(
           era == "early"  ~ sample(2007:2015, nrow(.), replace = TRUE),
           era == "middle" ~ sample(2016:2019, nrow(.), replace = TRUE),
           era == "recent" ~ sample(2020:2023, nrow(.), replace = TRUE)
         ))

# Parameters for simulation
# DV ~ bernoulli(p)
# logit(p) = alpha + beta_E + beta_S + (beta_T + beta_Et + beta_St) * T2007
# Where T2007 is "years since 2007"

# betas (E, S, Et, St)
beta_E <- c(marine = 0, terrestrial = 0.2)
beta_S <- c(aspatial = 0, spatial = 0.4)
beta_ET <- c(marine = 0, terrestrial = 0.03)
beta_ST <- c(aspatial = 0, spatial = 0.06)
beta_T <- 0.17
alpha <- -4

# Run simulation
inv_logit <- \(x) exp(x) / (1 + exp(x))
biolog <- biolog_grid %>%
  mutate(
    T2007 = T - 2007,
    logit_p = alpha + beta_E[E] + beta_S[S] + (beta_T + beta_ET[E] + beta_ST[S]) * T2007,
    p = inv_logit(logit_p),
    F = rbinom(nrow(.), size = 1, prob = p)
  )

# Visualize data
e_palette <- c(marine = "navy", terrestrial = "goldenrod")

biolog_text <- filter(biolog, T == 2023) %>%
  distinct(E, S, p) %>%
  mutate(label = str_to_title(paste(E, S, sep = "\n")))

ggplot(biolog, aes(T, p, color = E)) +
  geom_jitter(aes(y = F, shape = S),
              width = 0.25, height = 0.05,
              alpha = 0.5) +
  geom_line(aes(linetype = S)) +
  geom_text(aes(label = label), biolog_text,
            x = 2023.1, hjust = 0) +
  scale_color_manual(values = e_palette) +
  scale_x_continuous(limits = c(NA, 2025),
                     breaks = seq(2010, 2020, by = 5)) +
  labs(y = expression(p(f))) +
  expand_limits(x = 2025) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank())

# Fit model
biolog_prior <- c(
  set_prior(prior = "normal(0, 2)", class = "Intercept"),
  set_prior(prior = "normal(0, 0.2)", coef = "Eterrestrial"),
  set_prior(prior = "normal(0, 0.2)", coef = "Eterrestrial:T2007"),
  set_prior(prior = "normal(0, 0.2)", coef = "Sspatial"),
  set_prior(prior = "normal(0, 0.2)", coef = "T2007:Sspatial"),
  set_prior(prior = "normal(0, 0.2)", coef = "T2007")
)

# Prior predictive check
n_draws <- 1000
biolog_prior_pred <- tibble(
  .iter = seq(n_draws),
  alpha = rnorm(n_draws, 0, 2),
  b_Et = rnorm(n_draws, 0, 0.2),
  b_EtT = rnorm(n_draws, 0, 0.2),
  b_Ss = rnorm(n_draws, 0, 0.2),
  b_SsT = rnorm(n_draws, 0, 0.2),
  b_T = rnorm(n_draws, 0, 0.2)
) %>%
  cross_join(biolog) %>%
  mutate(b_E = ifelse(E == "terrestrial", b_Et, 0),
         b_S = ifelse(S == "spatial", b_Ss, 0),
         b_ET = ifelse(E == "terrestrial", b_EtT, 0),
         b_ST = ifelse(S == "spatial", b_SsT, 0),
         logit_p = alpha + b_E + b_S + (b_ET + b_ST + b_T) * T2007,
         p = inv_logit(logit_p))

# Linear predictor (20 random iterations)
biolog_prior_pred %>%
  filter(.iter %in% sample(.iter, 20)) %>%
  ggplot(aes(T, logit_p)) +
  geom_line(aes(color = E,
                linetype = S,
                group = interaction(.iter, E, S))) +
  theme_classic()
# Probability predictor (20 random iterations)
biolog_prior_pred %>%
  filter(.iter %in% sample(.iter, 20)) %>%
  ggplot(aes(T, p)) +
  geom_line(aes(color = E,
                linetype = S,
                group = interaction(.iter, E, S))) +
  theme_classic()

# Fit model
biolog_mod <- brm(
  F ~ E * T2007 + S * T2007,
  data = biolog,
  family = bernoulli(link = "logit"),
  prior = biolog_prior,
  chains = 4,
  iter = 10000,
  seed = 6789
)

# Posterior predictions: linear predictor
biolog_linpred <- linpred_draws(biolog_mod, biolog)
biolog_linpred %>%
  group_by(T, E, S) %>%
  summarize(logit_p_mean = mean(.linpred),
            logit_p_lwr = quantile(.linpred, 0.25),
            logit_p_upr = quantile(.linpred, 0.75),
            .groups = "drop") %>%
  ggplot(aes(T, logit_p_mean)) +
  geom_ribbon(aes(ymin = logit_p_lwr, ymax = logit_p_upr,
                  fill = E, linetype = S),
              alpha = 0.5) +
  geom_line(aes(y = logit_p_mean, color = E, linetype = S),
            linewidth = 1.2) +
  geom_line(aes(y = logit_p, color = E, linetype = S),
            biolog,
            linewidth = 3) +
  scale_color_manual(values = e_palette) +
  scale_fill_manual(values = e_palette) +
  scale_x_continuous(breaks = seq(2010, 2020, by = 5)) +
  labs(y = expression(logit(p(f)))) +
  theme_classic() +
  theme(axis.title.x = element_blank())

# Posterior predictions: probabilities
biolog_prob <- linpred_draws(biolog_mod, biolog, transform = TRUE)
biolog_prob %>%
  group_by(T, E, S) %>%
  summarize(p_mean = mean(.linpred),
            p_lwr = quantile(.linpred, 0.25),
            p_upr = quantile(.linpred, 0.75),
            .groups = "drop") %>%
  ggplot(aes(T, p_mean)) +
  geom_ribbon(aes(ymin = p_lwr, ymax = p_upr,
                  fill = E, linetype = S),
              alpha = 0.5) +
  geom_line(aes(y = p_mean, color = E, linetype = S),
            linewidth = 1.2) +
  geom_line(aes(y = p, color = E, linetype = S),
            biolog,
            linewidth = 3) +
  scale_color_manual(values = e_palette) +
  scale_fill_manual(values = e_palette) +
  scale_x_continuous(breaks = seq(2010, 2020, by = 5)) +
  labs(y = expression(p(f))) +
  theme_classic() +
  theme(axis.title.x = element_blank())

biolog_draws <- as_draws_df(biolog_mod) %>%
  as_tibble()
param_names <- c("b_Intercept", "b_Eterrestrial", "b_T2007", "b_Sspatial",
                 "b_Eterrestrial:T2007", "b_T2007:Sspatial")
param_actual <- c(alpha, beta_E["terrestrial"], beta_T, beta_S["spatial"],
                  beta_ET["terrestrial"], beta_ST["spatial"])
plot_grid(
  plotlist = map2(
    param_names,
    param_actual,
    \(coef, actual) {
      ggplot(biolog_draws, aes(.data[[coef]])) +
        geom_density(fill = "cornflowerblue", color = NA) +
        geom_density(fill = NA, color = "navy", linewidth = 1.5) +
        geom_vline(xintercept = actual, color = "firebrick", linewidth = 1.5) +
        theme_classic()
    }),
  ncol = 2)
