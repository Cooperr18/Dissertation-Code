# USING THE TEST FOR UNBIAS-------------------------------------------------

library(dplyr)
library(signatselect)
library(ggplot2)

# install.packages("pak")
pak::pkg_install("benmarwick/signatselect")

unbiased_mutation <- function(N, mu, p_0, t_max, r_max) {
  # Create the output tibble
  output <- tibble(generation = rep(1:t_max, r_max), 
                   p = as.numeric(rep(NA, t_max * r_max)), 
                   run = as.factor(rep(1:r_max, each = t_max))) 
  
  for (r in 1:r_max) { # Controls the number of runs in the experiment
    population <- tibble(trait = sample(c("A", "B"), N, replace = TRUE, 
                                        prob = c(p_0, 1 - p_0)))
    # Add first generation's p for run r
    output[output$generation == 1 & output$run == r, ]$p <- 
      sum(population$trait == "A") / N 
    for (t in 2:t_max) {
      # Copy individuals to previous_population tibble
      previous_population <- population 
      
      # Determine the probability of 'mutant' individuals
      mutate <- sample(c(TRUE, FALSE), N, prob = c(mu, 1 - mu), replace = TRUE) 
      
      # If there are 'mutants' from A to B
      if (nrow(population[mutate & previous_population$trait == "A", ]) > 0) { 
        # Then flip them to B
        population[mutate & previous_population$trait == "A", ]$trait <- "B" 
      }
      
      # If there are 'mutants' from B to A
      if (nrow(population[mutate & previous_population$trait == "B", ]) > 0) { 
        # Then flip them to A
        population[mutate & previous_population$trait == "B", ]$trait <- "A" 
      }
      
      # Get p and put it into output slot for this generation t and run r
      output[output$generation == t & output$run == r, ]$p <- 
        sum(population$trait == "A") / N 
    }
  }
  # Export data from function
  output 
}

# Parameters for simulation
N <- 10000     # Population size
mu <- 0.05    # Mutation rate
p_0 <- 0.5    # Initial allele frequency
t_max <- 150  # Number of generations
r_max <- 1   # Number of replicate runs (use 1 to match `fit()` input format)


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
