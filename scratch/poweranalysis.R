library(brms)
library(cowplot)
library(tidybayes)
library(tidyverse)
set.seed(1234)

inv_logit <- function(x) exp(x) / (1 + exp(x))
e_palette <- c(marine = "navy", terrestrial = "goldenrod")


# Params ------------------------------------------------------------------

alpha <- -4
beta_E <- c(marine = 0, terrestrial = 0.5)
beta_S <- c(aspatial = 0, spatial = 1.5)
beta_T = 0.2


# Simulation --------------------------------------------------------------

sim_analysis <- function(nstar) {
  biolog <- expand_grid(
    era = c("early", "middle", "recent"),
    E = c("marine", "terrestrial"),
    S = c("aspatial", "spatial")
  ) %>%
    slice(rep(seq(nrow(.)), each = nstar)) %>%
    mutate(
      T = case_when(
        era == "early"  ~ sample(2007:2015, nrow(.), replace = TRUE),
        era == "middle" ~ sample(2016:2019, nrow(.), replace = TRUE),
        era == "recent" ~ sample(2020:2023, nrow(.), replace = TRUE)
      ),
      T2007 = T - 2007,
      logit_p = alpha + beta_E[E] + beta_S[S] + beta_T * T2007,
      p = inv_logit(logit_p),
      DV = rbinom(nrow(.), 1, p)
    )
  biolog_prior <- c(
    set_prior(prior = "normal(0, 4)", class = "Intercept"),
    set_prior(prior = "normal(0, 2)", coef = "Eterrestrial"),
    set_prior(prior = "normal(0, 2)", coef = "Sspatial"),
    set_prior(prior = "normal(0, 0.5)", coef = "T2007")
  )
  biolog_mod <- brm(
    DV ~ E + S + T2007,
    data = biolog,
    family = bernoulli(link = "logit"),
    prior = biolog_prior,
    chains = 4,
    cores = 4,
    iter = 50000,
    refresh = 0
  )
  biolog_actual <- expand_grid(T2007 = (2007:2023) - 2007,
                               E = c("marine", "terrestrial"),
                               S = c("aspatial", "spatial")) %>%
    mutate(logit_p = alpha + beta_E[E] + beta_S[S] + beta_T * T2007)
  biolog_draws <- biolog_mod %>%
    as_draws_df() %>%
    as_tibble() %>%
    cross_join(biolog_actual) %>%
    mutate(linpred_post = b_Intercept +
             b_Eterrestrial * (E == "terrestrial") +
             b_Sspatial * (S == "spatial") +
             b_T2007 * T2007)
  rmse <- with(biolog_draws, sqrt(mean((linpred_post - logit_p)^2)))
  viz <- biolog_draws %>%
    group_by(T2007, E, S) %>%
    summarize(
      across(
        linpred_post,
        list(mean = mean,
             lwr = \(x) HDInterval::hdi(linpred_post, 0.95)[1],
             upr = \(x) HDInterval::hdi(linpred_post, 0.95)[2])),
      .groups = "drop") %>%
    ggplot(aes(T2007, color = E, fill = E, linetype = S)) +
    geom_ribbon(aes(ymin = linpred_post_lwr, ymax = linpred_post_upr),
                color = NA, alpha = 0.5) +
    geom_line(aes(y = linpred_post_mean)) +
    scale_color_manual(values = e_palette) +
    scale_fill_manual(values = e_palette) +
    geom_line(aes(y = logit_p), biolog_actual, linewidth = 1.5) +
    labs(title = sprintf("N* = %d", nstar),
         caption = sprintf("RMS error = %0.3f", rmse)) +
    theme_classic()
  list(rmse = rmse,
       viz = viz)
}


# Power analysis ----------------------------------------------------------

nstar <- seq(10, 100, by = 10)
niter <- 20
rmse <- vector(mode = "list", length = length(nstar))
names(rmse) <- nstar
viz <- vector(mode = "list", length = length(nstar))
names(viz) <- nstar
for (n in nstar) {
  rmse[[as.character(n)]] <- vector(mode = "double", length = niter)
  viz[[as.character(n)]] <- vector(mode = "list", length = niter)
  for (i in seq(niter)) {
    result <- sim_analysis(n)
    rmse[[as.character(n)]][i] <- result$rmse
    ggsave(sprintf("scratch/poweranalysisfigs/nstar=%d_i=%d.png", n, i),
           result$viz)
    rm(result)
  }
}


