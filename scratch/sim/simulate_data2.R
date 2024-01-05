library(tidyverse)

# Simulate data -----------------------------------------------------------
# We are trying to estimate a latent variable, "openness", from multiple
# indicators (F1, F4, ..., Complete) and predictors (year, habitat, and sensor
# type).

# First simulate 200 papers
n_papers <- 200
set.seed(1234)
indicators <- c("F1", "F4", "A1.1", "I1", "R1.1", "Complete")
openbiolog <- expand_grid(
  year = 2007:2022,
  habitat = c("marine", "terr"),
  sensor = c("spatial", "aspatial")
) %>%
  sample_n(n_papers, replace = TRUE) %>%
  mutate(id = factor(seq(n_papers)),
         year2007 = year - 2007)

# Simulate the "openness" of each paper. The distribution of openness should
# have a long tail, so we'll use a lognormal
X <- model.matrix(~ year2007 + year2007:habitat + year2007:sensor,
                  data = openbiolog)

as_vert <- function(x) t(t(x))
betas <- as_vert(c(
  1,     # intercept
  0.1,   # year2007
  0.2,  # year2007:habitatterr
  0.1    # year2007:sensorspatial
))
open_logsd <- 0.25

openbiolog <- openbiolog %>%
  mutate(open_eta = drop(X %*% betas),
         open_mu = log(open_eta),
         openness = rlnorm(n_papers,
                           meanlog = open_mu,
                           sdlog = open_logsd))

# Openness increases over time
# Marine (blue) < terrestrial (orange)
ggplot(openbiolog, aes(year2007, openness)) +
  geom_point(aes(shape = sensor, color = habitat)) +
  geom_smooth(method = "loess",
              se = FALSE, formula = y ~ x,
              color = "grey75") +
  scale_color_manual(values = c("cornflowerblue", "tan1")) +
  expand_limits(y = 0) +
  theme_classic()

# Openness is our latent variable. Now simulate indicator variables.
# Indicators are correlated within papers (e.g., if a paper satisfies F1 then
# it's more likely F4 is satisfied, too) but with different expected values
# (e.g., F1 is on average more likely to be satisfied than F4). Also, indicators
# are more likely to be satisfied for greater "openness" values.
# To simulate these properties, first generate multivariate normal values for
# indicators by paper. These values are correlated, but unrelated to openness.
# Then generate random normal values with a mean equal to "openness" and sd=1.
# If the random openness number exceeds the multivariate normal number, then the
# indicator is satisfied.
# I fully acknowledge this can't be the right way to do this!!

# Means and variance-covariance matrix of indicators
indicator_mu <- c(
  F1 = 3,
  F4 = 5,
  A1.1 = 5,
  I1 = 5,
  R1.1 = 5,
  Complete = 5
)
indicator_sigma <- toeplitz(c(1.0, 0.5, 0.5, 0.5, 0.5, 0.2))
# Multivariate normal openness thresholds (independent of openness)
indicator_thr <- MASS::mvrnorm(n_papers, indicator_mu, indicator_sigma)
# Normal indicator values (dependent on openness)
indicator_val <- t(sapply(openbiolog$openness, \(x) rnorm(6, mean = x, sd = 1)))

# Simulate indicator values
openbiolog[paste(indicators, "thr", sep = "_")] <- indicator_thr
openbiolog[paste(indicators, "val", sep = "_")] <- indicator_val
openbiolog[indicators] <- indicator_val > indicator_thr

openbiolog %>%
  pivot_longer(all_of(indicators),
               names_to = "indicator",
               values_to = "indicator_value") %>%
  mutate(indicator_value = as.integer(indicator_value)) %>%
  ggplot(aes(year, indicator_value, color = habitat)) +
  geom_jitter(width = 0.1, height = 0.1,
              alpha = 0.5) +
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              formula = y ~ x) +
  facet_wrap(~ indicator) +
  scale_color_manual(values = c("cornflowerblue", "tan1")) +
  theme_classic()

cor(openbiolog[indicators])
