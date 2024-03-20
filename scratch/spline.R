library(brms)
library(cowplot)
library(splines)
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
# time b-spline
T_theta <- cumsum(c(0, 0.6, 1.1, 0.7, 0.4, 0.2))
n_T_knots <- 3
T_knots <- seq(2007, 2023, length.out = n_T_knots + 2)[-c(1, n_T_knots + 2)]
T_basis <- bs(biolog_grid$T, knots = T_knots, degree = 3)

# betas (E, S)
beta_E <- c(marine = 0, terrestrial = 0.5)
beta_S <- c(aspatial = 0, spatial = 1)
alpha <- -4

# Run simulation
inv_logit <- \(x) exp(x) / (1 + exp(x))
biolog <- biolog_grid %>%
  mutate(bsT = as.vector(T_basis %*% T_theta),
         logit_p = alpha + beta_E[E] + beta_S[S] + bsT,
         p = inv_logit(logit_p),
         F = rbinom(nrow(.), size = 1, prob = p))

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
  set_prior(prior = "normal(0, 10)", class = "Intercept"),
  set_prior(prior = "normal(0, 10)", coef = "Eterrestrial"),
  set_prior(prior = "normal(0, 10)", coef = "Sspatial"),
  set_prior(prior = "normal(0, 10)", coef = "bsTdegreeEQ3dfEQ3P31"),
  set_prior(prior = "normal(0, 10)", coef = "bsTdegreeEQ3dfEQ3P32"),
  set_prior(prior = "normal(0, 10)", coef = "bsTdegreeEQ3dfEQ3P33"),
  set_prior(prior = "normal(0, 10)", coef = "bsTdegreeEQ3dfEQ3P34"),
  set_prior(prior = "normal(0, 10)", coef = "bsTdegreeEQ3dfEQ3P35"),
  set_prior(prior = "normal(0, 10)", coef = "bsTdegreeEQ3dfEQ3P36")
)

# Prior predictive check
n_draws <- 1000
biolog_prior_pred <- tibble(
  .iter = seq(n_draws),
  alpha = rnorm(n_draws, 0, 10),
  Eterrestrial = rnorm(n_draws, 0, 10),
  Sspatial = rnorm(n_draws, 0, 10),
  bsT_theta = matrix(rnorm(n_draws * ncol(T_basis), 0, 100), nrow = n_draws)
) %>%
  cross_join(biolog_grid) %>%
  group_by(.iter) %>%
  mutate(bsT = as.vector(T_basis %*% bsT_theta[1, ]),
         logit_p = alpha +
           ifelse(E == "terrestrial", Eterrestrial, 0) +
           ifelse(S == "spatial", Sspatial, 0) +
           bsT,
         p = inv_logit(logit_p)) %>%
  ungroup()

# Linear predictor (10 random iterations)
biolog_prior_pred %>%
  filter(.iter %in% sample(.iter, 10)) %>%
  ggplot(aes(T, logit_p)) +
  geom_line(aes(color = E,
                linetype = S,
                group = interaction(.iter, E, S))) +
  theme_classic()
# Probability predictor (10 random iterations)
biolog_prior_pred %>%
  filter(.iter %in% sample(.iter, 10)) %>%
  ggplot(aes(T, p)) +
  geom_line(aes(color = E,
                linetype = S,
                group = interaction(.iter, E, S))) +
  theme_classic()

# Fit model
biolog_mod <- brm(
  F ~ E + S + bs(T, degree = 3, df = 3 + 3), # 3 knots + 3rd degree
  data = biolog,
  family = bernoulli(link = "logit"),
  prior = biolog_prior,
  chains = 4,
  iter = 5000,
  seed = 6789
)

# Posterior predictions: linear predictor
biolog_linpred <- tidybayes::linpred_draws(
  biolog_mod,
  biolog
)
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

biolog_prob <- tidybayes::linpred_draws(
  biolog_mod,
  biolog,
  transform = TRUE
)
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


# 2 knots, 2nd order polynomial
# Fit model
biolog_prior <- c(
  set_prior(prior = "normal(0, 10)", class = "Intercept"),
  set_prior(prior = "normal(0, 10)", coef = "Eterrestrial"),
  set_prior(prior = "normal(0, 10)", coef = "Sspatial"),
  set_prior(prior = "normal(0, 10)", coef = "bsTdegreeEQ2dfEQ2P2interceptEQFALSE1"),
  set_prior(prior = "normal(0, 10)", coef = "bsTdegreeEQ2dfEQ2P2interceptEQFALSE2"),
  set_prior(prior = "normal(0, 10)", coef = "bsTdegreeEQ2dfEQ2P2interceptEQFALSE3"),
  set_prior(prior = "normal(0, 10)", coef = "bsTdegreeEQ2dfEQ2P2interceptEQFALSE4")
)
biolog_mod <- brm(
  F ~ E + S + bs(T, degree = 2, df = 2 + 2, intercept = FALSE), # 2 knots + 22nd degree
  data = biolog,
  family = bernoulli(link = "logit"),
  prior = biolog_prior,
  chains = 4,
  iter = 5000,
  seed = 6789
)

# Posterior predictions: linear predictor
biolog_linpred <- tidybayes::linpred_draws(
  biolog_mod,
  biolog
)
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

biolog_prob <- tidybayes::linpred_draws(
  biolog_mod,
  biolog,
  transform = TRUE
)
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

biolog_draws <- as_draws_df(biolog_mod)

mod_basis <- bs(biolog$T, degree = 2, df = 2 + 2, intercept = FALSE)
plot_grid(
  ggplot(biolog_draws, aes(b_Intercept)) +
    geom_density(fill = "cornflowerblue", color = NA) +
    geom_density(fill = NA, color = "navy", linewidth = 1.5) +
    geom_vline(xintercept = alpha, color = "firebrick", linewidth = 1.5) +
    theme_classic(),
  ggplot(biolog_draws, aes(b_Eterrestrial)) +
    geom_density(fill = "cornflowerblue", color = NA) +
    geom_density(fill = NA, color = "navy", linewidth = 1.5) +
    geom_vline(xintercept = alpha, color = "firebrick", linewidth = 1.5) +
    theme_classic(),
  ggplot(biolog_draws, aes(b_Sspatial)) +
    geom_density(fill = "cornflowerblue", color = NA) +
    geom_density(fill = NA, color = "navy", linewidth = 1.5) +
    geom_vline(xintercept = alpha, color = "firebrick", linewidth = 1.5) +
    theme_classic(),
  biolog_draws %>%
    as_tibble() %>%
    pivot_longer(starts_with("b_bsTdegree"),
                 names_to = "i",
                 values_to = "theta") %>%
    group_by(.chain, .iteration, .draw) %>%
    reframe(T = biolog$T,
            b_T = as.vector(mod_basis %*% theta)) %>%
    group_by(T) %>%
    summarize(b_T_mean = mean(b_T),
              b_T_lwr = quantile(b_T, 0.25),
              b_T_upr = quantile(b_T, 0.75)) %>%
    ggplot(aes(T)) +
    geom_ribbon(aes(ymin = b_T_lwr, ymax = b_T_upr),
                fill = "cornflowerblue") +
    geom_line(aes(y = b_T_mean),
              color = "navy",
              linewidth = 1.5) +
    geom_line(aes(y = bsT),
              biolog,
              color = "firebrick",
              linewidth = 1.5) +
    theme_classic(),
  nrow = 2
)

