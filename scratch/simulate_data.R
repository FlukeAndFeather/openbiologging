library(SimCorMultRes)
library(tidyverse)

# Simulated data to use as raw material for a discussion about methods
# Data will have the following properties
# Six response variables, call them f1, f4, a1.1, i1, r1.1, complete
# Six predictor variables: year, taxa, habitat, sensor type, biologging or other
# 250 papers evenly distributed across predictor variables
# All papers include biologging data, 50% of papers include other data
# Response variables are correlated. p(f1) is greatest; p(f4), p(a1.1), and
# p(r1.1) are lower and correlated with each other; p(i1) is lowest and
# correlated with p(f4); p(complete) is correlated with p(i1) and p(f4)

set.seed(789)
n_papers <- 250
unique_taxa <- c("mammal", "bird", "fish", "other")
unique_habitat <- c("marine", "terrestrial", "aquatic")

# Generating correlated random binomial variables
# The correlation matrix between the response variables are in a csv file
# Given six response variables, 2^6=64 potential combinations of 0's and 1's
# If we encode each paper's responses as 000000, 000001, 000010, 000011, etc.
# then we can randomly sample the numbers 1-64 with different probabilities for
# each number based on the correlations.
corr_tbl <- read_csv("scratch/sim/cormat.csv")
corr_matrix <- as.matrix(corr_tbl[, -1])
rownames(corr_matrix) <- corr_tbl[[1]]
corr_matrix

# hypothetical marginal probabilities
p <- c(f1 = 0.5,
       f4 = 0.3,
       a1.1 = 0.3,
       i1 = 0.15,
       r1.1 = 0.3,
       complete = 0.15)
p

p.joint <- mipfp::ObtainMultBinaryDist(corr = corr_matrix, marg.probs = p)

y.sim <- mipfp::RMultBinary(n = n_papers, mult.bin.dist = p.joint)
y.sim$binary.sequences

sim_data_biologging <- tibble(
  paperid = seq(n_papers),
  year = sample(2007:2022, n_papers, replace = TRUE),
  taxa = sample(unique_taxa, n_papers, replace = TRUE),
  habitat = sample(unique_habitat, n_papers, replace = TRUE),
  sensor = sample(c("spatial", "aspatial"), n_papers, replace = TRUE)
) %>%
  cbind(y.sim$binary.sequences) %>%
  pivot_longer(f1:complete, names_to = "rubric", values_to = "outcome")

fair_model <- geepack::geeglm(outcome ~ rubric + year + taxa + habitat + sensor,
                              data = sim_data_biologging,
                              id = paperid,
                              family = binomial)

summary(fair_model)


sim_biologging_data <- simulated_predictors %>%
  cbind(foo$ydata) %>%
  rename()






set.seed(123)
simulated_predictors <- tibble(
  paperid = seq(n_papers),
  year = sample(2007:2022, n_papers, replace = TRUE),
  taxa = sample(unique_taxa, n_papers, replace = TRUE),
  habitat = sample(unique_habitat, n_papers, replace = TRUE),
  sensor = sample(c("spatial", "aspatial"), n_papers, replace = TRUE)
)
sample_size <- 250
cluster_size <- 6
beta_intercepts <- runif(cluster_size)
beta_coefficients <- runif(
  1 + # year
    n_distinct(simulated_predictors$taxa) - 1 +
    n_distinct(simulated_predictors$habitat) - 1 +
    n_distinct(simulated_predictors$sensor) - 1
)
latent_correlation_matrix <- toeplitz(c(1, rep(0.75, cluster_size - 1)))
simulated_dataset <- rbin(
  clsize = cluster_size,
  intercepts = beta_intercepts,
  betas = beta_coefficients,
  xformula = ~ year + taxa + habitat + sensor,
  xdata = select(simulated_predictors, -paperid),
  cor.matrix = latent_correlation_matrix,
  link = "logit"
)

logit <- function(p) log(p / (1 - p))
logistic <- function(alpha) (tanh(alpha / 2) + 1) / 2

f1_year <- function(year, taxa, habitat, sensor) {
  beta_f1 <- c(
    intercept = logit(0.15),
    year = 0.1,
    taxa_mammal = 0,
    taxa_bird = 0.1,
    taxa_fish = -0.05,
    taxa_other = 0,
    habitat_marine = 0,
    habitat_terrestrial = 0.15,
    habitat_aquatic = -0.1,
    sensor_spatial = 0,
    sensor_aspatial = -0.1
  )
  X <- rbind(rep(1, length(year)),
             year,
             taxa == "mammmal",
             taxa == "bird",
             taxa == "fish",
             taxa == "other",
             habitat == "marine",
             habitat == "terrestrial",
             habitat == "aquatic",
             sensor == "spatial",
             sensor == "aspatial")
  rbinom(length(year),
         size = 1, # Bernoulli
         prob = logistic(t(beta_f1) %*% X))
}
simulated_dataset <- tibble(
  paperid = seq(n_papers),
  year = sample(2007:2022, n_papers, replace = TRUE),
  taxa = sample(unique_taxa, n_papers, replace = TRUE),
  habitat = sample(unique_habitat, n_papers, replace = TRUE),
  sensor = sample(c("spatial", "aspatial"), n_papers, replace = TRUE)
) %>%
  mutate(year_2007 = year - 2007,
         f1 = f1_year(year_2007, taxa, habitat, sensor))








