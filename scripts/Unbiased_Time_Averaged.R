# UNBIASED TRANSMISSION "TIME AVERAGED" VERSION ----------------------------------

library(signatselect)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(purrr)

set.seed(1234)

# Parameters --------
N <- 100 # Number of individuals
mu <- 0.2 # innovation rate (numeric, between 0 and 1, inclusive)
burnin <- 1000 # number of initial steps (iterations) discarded
timesteps <- 1000 # actual number of time steps or "generations" after the burn-in
p_value_lvl <- 0.05 # Significance level
n_runs <- 200 # number of test runs
time_window <- 20 # size of averaging windows

neutral_counts_per_run_ta <- numeric(n_runs) # empty vector for counting neutral variants
accuracy_ta <- numeric(n_runs) # empty vector for accuracy tracking each run

# Pipeline ---------

for (run in 1:n_runs) {
  ini <- 1:N # Initial cultural variants
  traitmatrix <- matrix(NA,nrow=timesteps,ncol=N) # Each row is a time step, and each column an individual
  pop <- ini # initial population of variants to pop
  maxtrait <- N 
  
  # Burn-in stage
  for(i in 1:burnin) {
    # neutral transmission:
    pop <- sample(pop,replace=T)
    
    # Add innovations
    innovate <- which(runif(N)<mu) # which individual innovates (>0.01 = TRUE)
    if(length(innovate) > 0) {
      new_variants <- (maxtrait + 1):(maxtrait + length(innovate)) # add variants
      pop[innovate] <- new_variants
      maxtrait <- max(pop)
    }
  }
  
  # Observation period after equilibrium
  for (i in 1:timesteps) {
    pop <- sample(pop,replace=T)
    innovate <- which(runif(N)<mu) 
    if(length(innovate) > 0) {
      new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
      pop[innovate] <- new_variants
      maxtrait <- max(pop)
    }
    
    traitmatrix[i,] <- pop #record the variants of each individual
  }
  
  # Time averaging step -----------------
  averaged_rows <- floor(timesteps / time_window) # round down to nearest whole number
  # list of time-averaged samples (each is a vector of N * time_window variants)
  averaged_samples <- vector("list", averaged_rows)
  for (j in 1:averaged_rows) {
    start <- (j - 1) * time_window + 1
    end <- j * time_window
    # Flatten all traits in the time window into one vector
    averaged_samples[[j]] <- as.vector(traitmatrix[start:end, ])
  }
  
  unique_variants <- sort(unique(unlist(averaged_samples))) # store unique variants across all bins
  
  # Trait matrix to frequency matrix (row = bins, col = variants)
  freq_mat <- t(sapply(averaged_samples, function(traits) {
    tab <- table(factor(traits, levels = unique_variants))
    as.numeric(tab) / length(traits)  # Proportions relative to N * time_window
  }))
  colnames(freq_mat) <- unique_variants
  
  # Prepare FIT input
  freq_long <- as.data.frame(freq_mat) %>%
    mutate(time = 1:nrow(.)) %>%
    pivot_longer(-time, names_to="variant", values_to="freq") %>% # long format
    filter(freq > 0) %>% # remove zeros
    mutate(variant = as.integer(variant))
  
  # Filter data to 3 time points
  freq_long_filtered <- freq_long %>%
    group_by(variant) %>%
    filter(n_distinct(time) >= 3) %>%
    ungroup()
  
  # FIT application and storing
  fit_results <- freq_long_filtered %>%
    group_split(variant) %>%
    map_dfr(~ {
      df <- as.data.frame(.x)  # .x is a dataframe
      
      # check if we have enough data points
      if (nrow(df) < 3) {
        return(data.frame(variant = df$variant[1], time_points = nrow(df), fit_p = NA, stringsAsFactors = FALSE))
      }
      
      # safely apply FIT
      res <- tryCatch(
        fit(time = df$time, v = df$freq),
        error = function(e) list(fit_p = NA)  # handle errors
      )
      
      data.frame(
        variant = df$variant[1],
        time_points = nrow(df),
        fit_p = res$fit_p,
        stringsAsFactors = FALSE
      )
    }) %>%
    mutate(
      sig = ifelse(fit_p > p_value_lvl, "neutral", "selection"), # neutral if >0.05
      sig = ifelse(is.na(fit_p), "NA", sig)  # missing values
    )
  
  # Store metrics
  total_variants <- nrow(fit_results)
  FPR <- sum(fit_results$sig == "selection") / total_variants  # False positives
  TNR <- sum(fit_results$sig == "neutral") / total_variants     # True negatives
  
  accuracy_ta[run] <- TNR  # track TNR across runs
}

# Check results
FPR
TNR
accuracy_ta[70] # we can check each run individually
overall_accuracy <- mean(accuracy_ta, na.rm = TRUE) # mean accuracy across runs
overall_accuracy

# Proportion of the runs have a 95% of detection
over95 <- sum(accuracy_ta >= 0.95)
high_accuracy_runs <- over95/n_runs*100
high_accuracy_runs

# Proportion of NA from the simulation
sumNA <- sum(fit_results$sig == "NA")
sumNEUTRAL <- sum(fit_results$sig == "neutral")
percentageNA <- sumNA/sumNEUTRAL * 100
percentageNA

# PLOTS --------------------------------------------

# Plot distribution of TNR, marking the 95% threshold
plot_neutral_ta <- ggplot(data.frame(TNR = accuracy_ta), aes(x = TNR)) +
  geom_histogram(binwidth = 0.002, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Neutral Detection Rate (NDR) Across Runs 'Time Averaged' Model", 
       subtitle = "Red line = expected TNR (1 - α)", 
       x = "True Neutral Rate", 
       y = "Frequency",
       caption = paste("Mean =", round(mean(accuracy_ta), 3), "|", 
                       "Runs ≥ 95% =", round(high_accuracy_runs, 1), "%", "|",
                       "Number of runs =", n_runs, "|",
                       "% NA =", round(percentageNA, 2), "%")) +
  theme_minimal()

grid.arrange(plot_neutral_snapshot,plot_neutral_ta,ncol=1)
plot_neutral_ta


