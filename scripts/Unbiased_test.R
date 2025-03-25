# UNBIASED TRANSMISSION ----------------------------------

remotes::install_github("benmarwick/signatselect", INSTALL_opts = "--no-multiarch --no-test-load")
library(signatselect)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# The parameters of the model are:
#   N: Population size (integer)
#   mu: innovation rate (numeric, between 0 and 1, inclusive)
#   burnin: number of initial steps (iterations) discarded
#   timesteps: actual number of time steps or "generations" after the burnin

N <- 100 # Number of traits
mu <- 0.01
burnin <- 100 # Get rid of arbitrary initial conditions
timesteps <- 100

# Significance level
p_value_lvl <- 0.05

# Run the test 100 times and store the results
n_runs <- 100
neutral_counts_per_run <- numeric(n_runs)

for (run in 1:n_runs) {
  ini <- 1:N # Initial cultural variants
  traitmatrix <- matrix(NA,nrow=timesteps,ncol=N) # Each row is a time step (100), and each column a trait (100)

  pop <- ini # initial population of variants to pop
  maxtrait <- N 

  # Loops before reaching equilibrium, determined by "burnin"
  for(i in 1:burnin) {
    # neutral transmission:
    pop <- sample(pop,replace=T)
  
    # Add innovations
    innovate <- which(runif(N)<mu) # which individual innovates (>0.01 = TRUE)
    if (length(innovate)>0) #if anyone innovates
    {
    number.new.variants <- length(innovate) # record the number
    new.variants <- maxtrait:(maxtrait+number.new.variants - 1) # create new variants
    pop[innovate] <- new.variants # assign the innovations to the innovators
    maxtrait <- max(pop) # update maxtrait matrix
    }
  }

  # And loops after reaching equilibrium
  for (i in 1:timesteps)
  {
    #transmission:
    pop <- sample(pop,replace=T)
    #innovation
    innovate <- which(runif(N)<mu) 
    if (length(innovate)>0) 
    {
      number.new.variants <- length(innovate) #how many innovators?
      new.variants <- maxtrait:(maxtrait+number.new.variants - 1) #create new variants
      pop[innovate] <- new.variants #assign the innovations to the innovators
      maxtrait <- max(pop) #update maxtrait 
    }
  
    traitmatrix[i,] <- pop #record the variants of each individual
  }

  unique_values <- sort(unique(as.vector(traitmatrix))) #all unique variants observed

  freq_mat <- t(apply(traitmatrix, 1, function(row) {
    tab <- table(factor(row, levels = unique_values))  # Count occurrences, ensuring all values are represented
    as.numeric(tab)
  }))

  colnames(freq_mat) <- unique_values # Assign column names for readability
  # Prepare input for fit() function
  freq.df <- as.data.frame(freq_mat)
  freq.df$time <- 1:timesteps # add time column

  #Convert to long format for fit()
  freq_long <- freq.df %>%
    pivot_longer(cols = -time, names_to = "variant", values_to = "freq") %>%
    mutate(variant = as.integer(variant))  # Variant as integers
  
  # Compute relative freq
  freq_long <- freq_long %>%
    group_by(time) %>%
    mutate(freq = freq / sum(freq)) %>%
    ungroup()
  
  # Safe version of fit, so that it doesn't stop, it just returns NA
  fit_safely <- safely(fit, 
                       otherwise = data.frame(fit_stat = NA, 
                                              fit_p = NA))
  
  # Filter variants with less than 3 time points
  list_of_dfs_three_or_more <- keep(list_of_dfs, ~sum(.x$freq > 0) >= 3)
  
  df_fit_test_results <- list_of_dfs_three_or_more %>%
    bind_rows(.id = "variant") %>%
    nest(-variant) %>%
    mutate(fit_test = map(data, ~fit_safely(time = .x$time, v = .x$freq))) %>%
    mutate(fit_p = map(fit_test, ~.x$result %>% bind_rows)) %>%
    unnest(fit_p) %>%
    mutate(sig = ifelse(fit_p > p_value_lvl, "neutral", "selection"))
  
  # Count neutral cases for this run
  neutral_counts_per_run[run] <- sum(df_fit_test_results$sig == "neutral", na.rm = TRUE)
}

# Final count of how often the FIT test detected neutrality
total_neutral_detections <- sum(neutral_counts_per_run)

cat("Total neutral detections out of", (n_runs * N), "tests:", total_neutral_detections, "\n")
cat("Proportion of neutral detections:", total_neutral_detections / (n_runs * N), "\n")


# --- PLOTS ---

# One single plot
ggplot(freq_long, aes(x = time, y = freq, group = variant, color = as.factor(variant))) +
  geom_line(alpha = 0.8) +
  labs(title = "Frequency Trajectories of All Variants", x = "Time", y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "none")

# Facet
ggplot(freq_long, aes(x = time, y = freq, group = variant)) +
  geom_line() +
  facet_wrap(~variant) +
  theme_minimal(base_size = 8) +
  labs(x = "Time", y = "Frequency") +
  theme(legend.position = "none")

ggplot(freq_long_sig, aes(time, freq, colour = sig, group = variant)) +
  geom_smooth(se = FALSE, method = "loess") +  # Smooth instead of raw lines
  facet_wrap(~variant, scales = "free_y") +
  theme_minimal()

                   
                   
                   
                   