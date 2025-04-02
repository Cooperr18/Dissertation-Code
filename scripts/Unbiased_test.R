# UNBIASED TRANSMISSION "SNAPSHOT" VERSION ----------------------------------

remotes::install_github("benmarwick/signatselect", INSTALL_opts = "--no-multiarch --no-test-load")
library(signatselect)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

N <- 80 # Number of traits
mu <- 0.01 # innovation rate (numeric, between 0 and 1, inclusive)
burnin <- 100 # number of initial steps (iterations) discarded
timesteps <- 100 # actual number of time steps or "generations" after the burn-in
p_value_lvl <- 0.05 # Significance level
n_runs <- 100 # number of test runs

neutral_counts_per_run_snapshot <- numeric(n_runs)

# Pipeline -----

for (run in 1:n_runs) {
  ini <- 1:N # Initial cultural variants
  traitmatrix <- matrix(NA,nrow=timesteps,ncol=N) # Each row is a time step (100), and each column a trait (100)
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
  
  unique_variants <- sort(unique(as.vector(traitmatrix)))

  # Trait matrix to frequency matrix
  freq_mat <- t(apply(traitmatrix, 1, function(row) {
    tab <- table(factor(row, levels = unique_variants))
    as.numeric(tab)/N # Convert to frequencies
  }))
  colnames(freq_mat) <- unique_variants

  # Prepare FIT input
  freq_long <- as.data.frame(freq_mat) %>%
    mutate(time = 1:timesteps) %>%
    pivot_longer(-time, names_to="variant", values_to="freq") %>% # long format
    filter(freq > 0) %>% # remove zeros
    mutate(variant = as.integer(variant))
  
  # Filter data to 3 time points
  freq_long_filtered <- freq_long %>%
    group_by(variant) %>%
    filter(n_distinct(time) >= 3) %>%
    ungroup()
  
  # Apply FIT
  fit_results <- freq_long_filtered %>%
    group_split(variant) %>%
    map_dfr(~ {
      df <- as.data.frame(.x)  # Ensure .x is a dataframe
      
      # Check if we have enough data points
      if (nrow(df) < 3) {
        return(data.frame(variant = df$variant[1], time_points = nrow(df), fit_p = NA, stringsAsFactors = FALSE))
      }
      
      # Safely apply the FIT test
      res <- tryCatch(
        fit(time = df$time, v = df$freq),
        error = function(e) list(fit_p = NA)  # Handle errors gracefully
      )
      
      data.frame(
        variant = df$variant[1],
        time_points = nrow(df),
        fit_p = res$fit_p,
        stringsAsFactors = FALSE
      )
    }) %>%
    mutate(
      sig = ifelse(fit_p > p_value_lvl, "neutral", "selection"),
      sig = ifelse(is.na(fit_p), "NA", sig)  # missing values
    )
  
  # Count neutral cases for a run
  neutral_counts_per_run_snapshot[run] <- sum(fit_results$sig == "neutral", na.rm = TRUE)
}

# Final count
total_neutral_detections_snapshot <- sum(neutral_counts_per_run_snapshot)
cat("Total neutral detections out of", (n_runs * N), "tests:", total_neutral_detections_snapshot, "\n")
cat("Proportion of neutral detections:", total_neutral_detections_snapshot / (n_runs * N), "\n")

# Visualization
freq_long_sig <- freq_long %>%
  left_join(fit_results %>% select(variant, sig), by = "variant")

unique(freq_long_filtered$variant)

# --- PLOTS ---

# One single plot
ggplot(freq_long_sig, aes(x = time, y = freq, color = sig, group = variant)) +
  geom_line() +
  scale_color_viridis_d(name = "", begin = 0.25, end = 0.75) +
  theme_minimal(base_size = 10) +
  ggtitle("FIT Classification of Neutral Transmission Data")

# Facet
ggplot(freq_long_filtered, aes(x = time, y = freq, group = variant)) +
  geom_line(stat = "identity") +  # lines are drawn for multiple points
  facet_wrap(~variant, scales = "free_y") +  # test is computed with total frequency, not relative
  theme_minimal(base_size = 8) +
  labs(x = "Time", y = "Frequency") +
  theme(legend.position = "none")





                   
                   
                   
                   