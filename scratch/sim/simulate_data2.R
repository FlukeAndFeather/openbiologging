library(ggeffects)
library(lme4)
library(tidyverse)

# Utility functions -------------------------------------------------------

logit <- function(p) log(p / (1 - p))
encode_dummy <- function(df, cols) {
  encode_one <- function(.df, col) {
    values <- sort(unique(df[[col]]))
    for (v in values[-1]) {
      .df[[paste(col, v, sep = "_")]] <- ifelse(df[[col]] == v, 1, 0)
    }
    .df
  }
  reduce(cols, encode_one, .init = df)
}
sim_ranef <- function(effect, effect_var) {
  n_effects <- length(unique(effect))
  values <- rnorm(n = n_effects, mean = 0, sd = sqrt(effect_var))
  names(values) <- unique(effect)
  values[effect]
}
predict_logodds <- function(int, coef, ...) {
  x <- do.call(rbind, list(...))
  as.vector(int + t(t(coef) %*% x))
}

# Simulation parameters ---------------------------------------------------

n_papers <- 100

set.seed(1234)

# "openness" intercept in log odds
intercept <- logit(0.1)

# fixed effect coefficients
beta <- c(
  year2007 = 0.1,
  habitat_terr = 0.2,
  taxa_fish = -0.3,
  taxa_mammal = -0.05,
  taxa_other = -0.15,
  sensor_spatial = 0.75
)

# random effects
id_var <- 1
fair_int_var <- 0.5
fair_year_var <- 1

set.seed(1234)
n_papers <- 100 # number of papers
n_scores <- 6   # number of measurements per paper

biologging_data <- expand_grid(
  year = 2007:2022,
  habitat = c("marine", "terr"),
  taxa = c("bird", "fish", "mammal", "other"),
  sensor = c("spatial", "aspatial")
) %>%
  sample_n(n_papers, replace = TRUE) %>%
  mutate(id = factor(seq(n_papers)),
         year2007 = year - 2007) %>%
  cross_join(tibble(attribute = c("F1", "F4", "A1.1", "I1", "R4", "Complete")))

X <- model.matrix(~ year2007 + habitat + taxa + sensor, data = biologging_data)
Z <- model.matrix(~ id + attribute, data = biologging_data)

betas <- c(
  logit(0.1), # Intercept
  0.1,        # year2007
  0.2,        # habitatterr
  -0.3,       # taxafish
  -0.05,      # taxamammal
  -0.15,      # taxaother
  0.75        # sensorspatial
)

# variance/covariance of random effects
# variance of random intercept for id
D_id <- 0.4
# variance/covariance of attributes
D_attr <- toeplitz(c(0.6, rep(0.3, 5)))

L <- chol(D_attr)
nvars <- dim(L)[1]
r = t(L) %*% matrix(rnorm(nvars), nrow = nvars, ncol = 1)

# Simulate data -----------------------------------------------------------

openbiologging <- expand_grid(
  year = 2007:2022,
  habitat = c("marine", "terr"),
  taxa = c("bird", "fish", "mammal", "other"),
  sensor = c("spatial", "aspatial"),
  fair_attr = c("F1", "F4", "A1.1", "I1", "R4", "Complete"),
  fair_val = 1
) %>%
  pivot_wider(names_from = fair_attr, values_from = fair_val) %>%
  sample_n(n_papers, replace = TRUE) %>%
  mutate(id = factor(seq(n_papers)),
         year2007 = year - 2007) %>%
  pivot_longer(F1:Complete, names_to = "fair", values_to = "drop") %>%
  select(-drop) %>%
  encode_dummy(c("habitat", "taxa", "sensor", "fair")) %>%
  mutate(fair_intercept = sim_ranef(fair, fair_int_var),
         fair_year_slope = sim_ranef(fair, fair_year_var),
         id_effect = sim_ranef(id, id_var),
         int = intercept + fair_intercept + id_effect,
         log_odds_fixed = predict_logodds(intercept, beta, year2007,
                                          habitat_terr, taxa_fish, taxa_mammal,
                                          taxa_other, sensor_spatial),
         log_odds_random = predict_logodds(fair_intercept + id_effect,
                                           beta, fair_F1, fair_),
         prop = plogis(log_odds),
         fair_val = rbinom(nrow(.), 1, prop))

# Visualize data ----------------------------------------------------------

plot_data <- function(predictor) {
  data_long <- openbiologging %>%
    count(.data[[predictor]], fair, fair_val)
  summary_geom <- if (predictor == "year") {
    geom_smooth(aes(weight = n),
                formula = y ~ x,
                method = glm,
                method.args = list(family = "binomial"))
  } else {
    geom_pointrange(aes(y = p, ymin = p - err, ymax = p + err),
                    data = data_long %>%
                      group_by(.data[[predictor]], fair) %>%
                      summarize(p = n[fair_val == 1] / sum(n),
                                err = sqrt(p * (1 - p) / sum(n)),
                                .groups = "drop"))
  }
  ggplot(data_long, aes(.data[[predictor]], fair_val)) +
    geom_point(aes(size = n), shape = 21) +
    summary_geom +
    facet_wrap(~ fair) +
    theme_classic()
}
# by year
plot_data("year")
plot_data("habitat")
plot_data("taxa")
plot_data("sensor")

# GLMM
# Doesn't converge w/ habitat + taxa
openbiologging_glmm <- glmer(
  cbind(fair_val, 1 - fair_val) ~ year2007 + habitat + taxa + sensor +
    (1 | id) + (year2007 | fair),
  data = openbiologging,
  family = binomial
)

# Plot predictions
openbiologging_pred <- expand_grid(
  year2007 = unique(openbiologging$year2007),
  habitat = unique(openbiologging$habitat),
  taxa = unique(openbiologging$taxa),
  sensor = unique(openbiologging$sensor),
  fair = unique(openbiologging$fair)
) %>%
  mutate(fair_logit = predict(openbiologging_glmm,
                              newdata = .,
                              type = "link",
                              re.form = ~ (1 | fair)),
         Year = year2007 + 2007)

# By habitat
openbiologging_pred %>%
  group_by(Year, habitat, fair) %>%
  summarize(fair_logit = mean(fair_logit), .groups = "drop") %>%
  mutate(fair_val = plogis(fair_logit)) %>%
  ggplot(aes(Year, fair_val)) +
  geom_line(aes(y = fair_val, color = habitat)) +
  facet_grid(cols = vars(fair)) +
  scale_y_continuous("p(FAIR)", limits = c(0, NA)) +
  theme_classic()

# By habitat
openbiologging_pred %>%
  group_by(Year, habitat, fair) %>%
  summarize(fair_logit = mean(fair_logit), .groups = "drop") %>%
  mutate(fair_val = plogis(fair_logit)) %>%
  ggplot(aes(Year, fair_val)) +
  geom_line(aes(y = fair_val, color = habitat)) +
  facet_grid(cols = vars(fair)) +
  scale_y_continuous("p(FAIR)", limits = c(0, NA)) +
  theme_classic()
