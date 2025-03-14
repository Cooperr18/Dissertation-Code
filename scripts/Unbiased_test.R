# UNBIASED TRANSMISSION PIPELINE ----------------------------------

library(dplyr)
library(signatselect)
library(ggplot2)
library(tibble)

pak::pkg_install("benmarwick/signatselect")

# The parameters of the simulation are:
#   N: Population size (integer)
#   mu: Mutation rate (numeric, between 0 and 1, inclusive)
#   p_0: Initial frequency of trait "A" (numeric, between 0 and 1, inclusive)
#   t_max: Number of generations (integer)
#   r_max: Number of replicate runs (integer)

unbiased_mutation <- function(N, mu, p_0, t_max, r_max) {
  # Validate p_0 and mu probabilities
  # If sampling returns negative or >1, we warn and provide a default value
  if(p_0 < 0 || p_0 > 1) {
    warning("p_0 must be between 0 and 1 (inclusive).")
    p_0 <- 0.5
  }
  if(mu < 0 || mu > 1) {
    warning("mu must be between 0 and 1 (inclusive).")
    mu <- 0.5
  }
  # Create the output tibble
  output <- tibble(generation = rep(1:t_max, r_max), 
                   p = as.numeric(rep(NA, t_max * r_max)), 
                   run = as.factor(rep(1:r_max, each = t_max))) 
  
  for (r in 1:r_max) { # Loop over runs
    # Initialize population of A and B based on p_0
    population <- tibble(trait = sample(c("A", "B"), N, replace = TRUE, 
                                        prob = c(p_0, 1 - p_0)))
    # Initial frequency of A
    output[output$generation == 1 & output$run == r, ]$p <- 
      sum(population$trait == "A") / N 
    
    for (t in 2:t_max) { # Loop over generations
      
      # Determine the probability of "mutant" individuals
      # If they mutate, TRUE is returned
      mutate <- sample(c(TRUE, FALSE), N, prob = c(mu, 1 - mu), replace = TRUE) 
      
      # Only select those traits that mutate
      population$trait[mutate] <- ifelse(population$trait[mutate] == "A", "B", "A")
      # mutants which A = TRUE, B will be returned
      # mutants which A = FALSE, A will be returned
      
      # Get p and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p <- 
        sum(population$trait == "A") / N 
    }
  }
  # Export data from function
  output 
}

# Parameters for simulation
N <- 10000     
mu <- 0.05   
p_0 <- 0.5    
t_max <- 150  
r_max <- 1  # 1 to match "fit()" input format


sim_data <- unbiased_mutation(N, mu, p_0, t_max, r_max) # Simulate neutrality
plot_multiple_runs(sim_data)

# Prepare input for fit() function
neutral_test <- tibble(
  time = sim_data$generation,
  freq = sim_data$p
)

# Apply the fit() function
fit_result <- fit(time = neutral_test$time, v = neutral_test$freq)
print(fit_result)




