library(brms)
library(cowplot)
library(tidyverse)
set.seed(1234)

biolog_grid <- expand_grid(
  era = c("early", "middle", "recent"),
  E = c("marine", "terrestrial"),
  S = c("aspatial", "spatial")
) %>%
  slice(rep(seq(nrow(.)), each = 100)) %>%
  mutate(T = case_when(
    era == "early"  ~ sample(2007:2015, nrow(.), replace = TRUE),
    era == "middle" ~ sample(2016:2019, nrow(.), replace = TRUE),
    era == "recent" ~ sample(2020:2023, nrow(.), replace = TRUE)
  ))

# Parameters for simulation
# delta
delta <- c(
  # dummy for 2007
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
beta_E <- c(marine = 0, terrestrial = 0.5)
beta_S <- c(aspatial = 0, spatial = 1)
beta_T <- 2.5
alpha <- -4

# Run simulation
inv_logit <- \(x) exp(x) / (1 + exp(x))
biolog <- biolog_grid %>%
  mutate(cumsum_delta = map_dbl(T, \(t) sum(delta[1:(t - 2007 + 1)])),
         logit_p = alpha + beta_E[E] + beta_S[S] + beta_T * cumsum_delta,
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

# Simulate priors
sim_priors <- function(biolog_prior, n = 1000) {
  prior_grid <- expand_grid(
    E = c("marine", "terrestrial"),
    S = c("aspatial", "spatial"),
    T = 2007:2023,
    F = c(0, 1)
  )

  prior_mod <- brm(F ~ E + S + mo(T),
                   data = prior_grid,
                   family = bernoulli(link = "logit"),
                   prior = biolog_prior,
                   chains = 1,
                   iter = 1000,
                   refresh = 0,
                   sample_prior = TRUE)
  biolog_prior_draws <- prior_draws(prior_mod) %>%
    as_draws_df()
  simo_moT1_cols <- str_subset(colnames(biolog_prior_draws), "simo_moT1")
  prior_grid %>%
    select(-F) %>%
    distinct() %>%
    cross_join(biolog_prior_draws) %>%
    mutate(cumsum_delta = apply(.[c("T", simo_moT1_cols)],
                                1,
                                \(row) {
                                  t <- row[1]
                                  delta <- c(0, row[-1])
                                  sum(delta[1:(t - 2007 + 1)])
                                }),
           logit_p = Intercept +
             ifelse(E == "terrestrial", b_Eterrestrial, 0) +
             ifelse(S == "spatial", b_Sspatial, 0) +
             bsp_moT * length(simo_moT1_cols) * cumsum_delta,
           p = inv_logit(logit_p))
}

summarize_priors <- function(prior_df) {
  prior_df %>%
    group_by(T, E, S) %>%
    summarize(logit_p_lwr = quantile(logit_p, 0.025),
              logit_p_upr = quantile(logit_p, 0.975),
              logit_p_mean = mean(logit_p),
              .groups = "drop")
}

# "Standard normal" priors
prior_stdnorm <- sim_priors(
  c(
    set_prior(prior = "normal(0, 1)", class = "Intercept"),
    set_prior(prior = "normal(0, 1)", coef = "Eterrestrial"),
    set_prior(prior = "normal(0, 1)", coef = "Sspatial"),
    set_prior(prior = "normal(0, 1)", coef = "moT"),
    set_prior(prior = "dirichlet(1)", class = "simo", coef = "moT1")
  )
)
ggplot(prior_stdnorm, aes(T, logit_p, group = .draw)) +
  geom_line(alpha = 0.1) +
  facet_grid(rows = vars(S), cols = vars(E)) +
  theme_classic()
ggplot(prior_stdnorm, aes(p, color = E, linetype = S)) +
  geom_density() +
  theme_classic()

# "Informative" priors
prior_inform <- sim_priors(
  c(
    set_prior(prior = "normal(0, 1)", class = "Intercept"),
    set_prior(prior = "normal(0, 0.5)", coef = "Eterrestrial"),
    set_prior(prior = "normal(0, 0.5)", coef = "Sspatial"),
    set_prior(prior = "normal(0, 0.1)", coef = "moT"),
    set_prior(prior = "dirichlet(1)", class = "simo", coef = "moT1")
  )
)
ggplot(prior_inform, aes(T, logit_p, group = .draw)) +
  geom_line(alpha = 0.1) +
  facet_grid(rows = vars(S), cols = vars(E)) +
  theme_classic()
ggplot(prior_inform, aes(p, color = E, linetype = S)) +
  geom_density() +
  theme_classic()

# "Uninformative" priors
prior_uninform <- sim_priors(
  c(
    set_prior(prior = "normal(0, 100)", class = "Intercept"),
    set_prior(prior = "normal(0, 100)", coef = "Eterrestrial"),
    set_prior(prior = "normal(0, 100)", coef = "Sspatial"),
    set_prior(prior = "normal(0, 100)", coef = "moT"),
    set_prior(prior = "dirichlet(1)", class = "simo", coef = "moT1")
  )
)
ggplot(prior_uninform, aes(T, logit_p, group = .draw)) +
  geom_line(alpha = 0.1) +
  facet_grid(rows = vars(S), cols = vars(E)) +
  theme_classic()
ggplot(prior_uninform, aes(p, color = E, linetype = S)) +
  geom_density() +
  theme_classic()
