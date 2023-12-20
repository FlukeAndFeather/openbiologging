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

# "openness" intercept in log odds
intercept <- logit(0.1)

# coefficients (betas)
beta <- c(
  year2007 = 0.15,
  habitat_terr = 1,
  taxa_fish = -1.5,
  taxa_mammal = -1,
  taxa_other = -2,
  sensor_spatial = 3
)

# random effect variances
fair_var <- 0.5
id_var <- 1

# Simulate data -----------------------------------------------------------

set.seed(1234)
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
  encode_dummy(c("habitat", "taxa", "sensor")) %>%
  mutate(fair_effect = sim_ranef(fair, fair_var),
         id_effect = sim_ranef(id, id_var),
         int = intercept + fair_effect + id_effect,
         log_odds = predict_logodds(int, beta, year2007, habitat_terr,
                                    taxa_fish, taxa_mammal, taxa_other,
                                    sensor_spatial),
         prop = plogis(log_odds))


# Visualize data ----------------------------------------------------------





openbiologging <- tibble(
  id = 1:n_papers,
  year = sample(2007:2022, n_papers, replace = TRUE),
  habitat = sample(c("Marine", "Terrestrial"), n_papers, replace = TRUE),
  taxa = sample(c("Mammal", "Bird", "Fish", "Other"), n_papers, replace = TRUE),
  year2007 = year - 2007,
  p = prob(year2007, habitat, taxa, psigma = 1),
  findable = rbinom(n_papers, size = 1, prob = p),
  accessible = rbinom(n_papers,
                      size = 1,
                      # If it's findable, it's more likely to be accessible
                      prob = ifelse(findable,
                                    logistic(logit(p) + rnorm(n_papers, -0.5, 2)),
                                    logistic(logit(p) + rnorm(n_papers, -1.5, 2)))),
  reusable = rbinom(n_papers,
                    size = 1,
                    # If it's findable, it's more likely to be reusable
                    prob = ifelse(findable,
                                  logistic(logit(p) + rnorm(n_papers, -1, 3)),
                                  logistic(logit(p) + rnorm(n_papers, -2, 3))))
)

plot_data <- function(predictor) {
  data_long <- openbiologging %>%
    pivot_longer(findable:reusable,
                 names_to = "attribute",
                 values_to = "value") %>%
    count(.data[[predictor]], attribute, value)
  summary_geom <- if (predictor == "year") {
    geom_smooth(aes(weight = n),
                formula = y ~ x,
                method = glm,
                method.args = list(family = "binomial"))
  } else {
    geom_pointrange(aes(y = p, ymin = p - err, ymax = p + err),
                    data = data_long %>%
                      group_by(.data[[predictor]], attribute) %>%
                      summarize(p = n[value == 1] / sum(n),
                                err = sqrt(p * (1 - p) / sum(n)),
                                .groups = "drop"))
  }
  ggplot(data_long, aes(.data[[predictor]], value)) +
    geom_point(aes(size = n), shape = 21) +
    summary_geom +
    facet_grid(cols = vars(attribute)) +
    theme_classic()
}
# by year
plot_data("year")
plot_data("habitat")
plot_data("taxa")

# GLMM
# Doesn't converge w/ habitat + taxa
biologging_long <- openbiologging %>%
  pivot_longer(findable:reusable, names_to = "attribute", values_to = "value")
openbiologging_glmm <- glmer(
  cbind(value, 1 - value) ~ attribute * year2007 + habitat + (1 | id),
  data = biologging_long,
  family = binomial
)

# Plot predictions
ggemmeans(openbiologging_glmm,
          terms = c("year2007", "habitat", "attribute")) %>%
  mutate(Year = x + 2007) %>%
  ggplot(aes(Year, predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2) +
  geom_line(aes(y = predicted, color = group)) +
  facet_grid(cols = vars(facet)) +
  scale_y_continuous("p(FAIR)", limits = c(0, NA)) +
  theme_classic()

summary(openbiologging_glmm)
