# HES PRESENTATION --------------------------------------

library(signatselect)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(purrr)

set.seed(1234)

# Parameters --------
N_pres <- 50 # Number of individuals
mu_pres <- 0.02 # innovation rate (numeric, between 0 and 1, inclusive)
burnin_pres <- 1000 # number of initial steps (iterations) discarded
timesteps_pres <- 1000 # actual number of time steps or "generations" after the burn-in
p_value_lvl_pres <- 0.05 # Significance level
n_runs_pres <- 500 # number of test runs

time_window_pres <- 20 # size of averaging windows

neutral_counts_per_run_snapshot <- numeric(n_runs_pres) # empty vector for counting neutral variants
accuracy_snapshot <- numeric(n_runs_pres) # empty vector for accuracy tracking each run

neutral_counts_per_run_ta <- numeric(n_runs_pres) # empty vector for counting neutral variants
accuracy_ta <- numeric(n_runs_pres) # empty vector for accuracy tracking each run



for (run in 1:n_runs_pres) {
  ini <- 1:N_pres  # Initial cultural variants
  traitmatrix <- matrix(NA, nrow = timesteps_pres, ncol = N_pres)  # Each row is a time step, each column an individual
  pop <- ini  # Initial population
  maxtrait <- N_pres
  
  # Burn-in stage
  for(i in 1:burnin_pres) {
    pop <- sample(pop, replace = TRUE)
    
    # Add innovations
    innovate <- which(runif(N_pres) < mu_pres)
    if(length(innovate) > 0) {
      new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
      pop[innovate] <- new_variants
      maxtrait <- max(pop)
    }
  }
  
  # Observation period
  for (i in 1:timesteps_pres) {
    pop <- sample(pop, replace = TRUE)
    innovate <- which(runif(N_pres) < mu_pres)
    if(length(innovate) > 0) {
      new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
      pop[innovate] <- new_variants
      maxtrait <- max(pop)
    }
    
    traitmatrix[i, ] <- pop
  }
  
  unique_variants <- sort(unique(as.vector(traitmatrix)))
  
  # Trait matrix to frequency matrix
  freq_mat <- t(apply(traitmatrix, 1, function(row) {
    tab <- table(factor(row, levels = unique_variants))
    as.numeric(tab) / N_pres
  }))
  colnames(freq_mat) <- unique_variants
  
  # Prepare FIT input
  freq_long <- as.data.frame(freq_mat) %>%
    mutate(time = 1:timesteps_pres) %>%
    pivot_longer(-time, names_to = "variant", values_to = "freq") %>%
    filter(freq > 0) %>%
    mutate(variant = as.integer(variant))
  
  # Filter to 3+ time points
  freq_long_filtered <- freq_long %>%
    group_by(variant) %>%
    filter(n_distinct(time) >= 3) %>%
    ungroup()
  
  # FIT application and storing
  fit_results <- freq_long_filtered %>%
    group_split(variant) %>%
    map_dfr(~ {
      df <- as.data.frame(.x)
      if (nrow(df) < 3) {
        return(data.frame(variant = df$variant[1], time_points = nrow(df), fit_p = NA, stringsAsFactors = FALSE))
      }
      
      res <- tryCatch(
        fit(time = df$time, v = df$freq),
        error = function(e) list(fit_p = NA)
      )
      
      data.frame(
        variant = df$variant[1],
        time_points = nrow(df),
        fit_p = res$fit_p,
        stringsAsFactors = FALSE
      )
    }) %>%
    mutate(
      sig = ifelse(fit_p > p_value_lvl_pres, "neutral", "selection"),
      sig = ifelse(is.na(fit_p), "NA", sig)
    )
  
  # Store metrics
  total_variants <- nrow(fit_results)
  FPR_pres <- sum(fit_results$sig == "selection") / total_variants
  NDR_pres <- sum(fit_results$sig == "neutral") / total_variants
  
  accuracy_snapshot[run] <- NDR_pres
}

# Check results
FPR_pres
NDR_pres
accuracy_snapshot[70] # we can check each run individually
overall_accuracy_pres <- mean(accuracy_snapshot, na.rm = TRUE) # mean accuracy across runs
overall_accuracy_pres

# Proportion of the runs have a 95% of detection
over95 <- sum(accuracy_snapshot >= 0.95)
high_accuracy_runs_pres <- over95/n_runs*100
high_accuracy_runs_pres

# Proportion of NA from the simulation
sumNA_pres <- sum(fit_results$sig == "NA")
sumNA_pres
proportionNA_pres <- sumNA_pres/length(fit_results$sig)*100
proportionNA_pres

#PLOT ----------------------------------------

# Prepare data
freq_long_sig <- freq_long %>%
  left_join(fit_results %>% select(variant, sig), by = "variant")

# Filter to keep only variants with max freq > 0.2
freq_plot_data <- freq_long_filtered %>%
  group_by(variant) %>%
  filter(max(freq) > 0.4) %>%
  ungroup()

mean_freq <- freq_plot_data %>%
  group_by(time) %>%
  summarise(freq = mean(freq), .groups = "drop") %>%
  mutate(variant = "mean")  # Dummy label for ggplot aesthetics

# Unique plot
library(RColorBrewer)
plot_variant_freq <- ggplot(freq_plot_data, aes(x = time, y = freq, group = variant, color = as.factor(variant))) +
  geom_line() +
  geom_line(data = mean_freq, aes(x = time, y = freq), 
            color = "black", size = 1.25, linetype = "solid") +
  scale_color_brewer(palette = "Dark2", name = "") +
  theme_minimal(base_size = 10)

plot_variant_freq

# Facet plot
plot_variant_freq_facet <- ggplot(freq_plot_data, aes(x = time, y = freq, group = variant)) +
  geom_line(stat = "identity") +
  facet_wrap(~variant) +
  theme_minimal(base_size = 8) +
  labs(x = "Time", y = "Frequency") +
  theme(legend.position = "none")

plot_variant_freq_facet


# Plot distribution of NDR, marking the 95% threshold
plot_neutral_snapshot <- ggplot(data.frame(NDR = accuracy_snapshot), aes(x = NDR)) +
  geom_histogram(binwidth = 0.005, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Neutral Detection Rate (NDR) Across Runs 'Snapshot' Model", 
       subtitle = "Red line = expected NDR (1 - α)", 
       x = "Neutral Detection Rate", 
       y = "Frequency",
       caption = paste("Average =", round(mean(accuracy_snapshot), 3), "|", 
                       "Runs ≥ 95% =", round(high_accuracy_runs, 1), "%", "|",
                       "Number of runs =", n_runs, "|",
                       "% NA =", round(percentageNA, 2), "%")) +
  theme_minimal()

plot_neutral_snapshot

# TIME AVERAGING -------------------------------------------------

for (run in 1:n_runs_pres) {
  ini <- 1:N_pres  # Initial cultural variants
  traitmatrix <- matrix(NA, nrow = timesteps_pres, ncol = N_pres)  # Each row = time step, each column = individual
  pop <- ini  # Initial population
  maxtrait <- N_pres
  
  # Burn-in stage
  for(i in 1:burnin_pres) {
    pop <- sample(pop, replace = TRUE)
    
    # Add innovations
    innovate <- which(runif(N_pres) < mu_pres)
    if(length(innovate) > 0) {
      new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
      pop[innovate] <- new_variants
      maxtrait <- max(pop)
    }
  }
  
  # Observation period after equilibrium
  for (i in 1:timesteps_pres) {
    pop <- sample(pop, replace = TRUE)
    innovate <- which(runif(N_pres) < mu_pres)
    if(length(innovate) > 0) {
      new_variants <- (maxtrait + 1):(maxtrait + length(innovate))
      pop[innovate] <- new_variants
      maxtrait <- max(pop)
    }
    
    traitmatrix[i, ] <- pop
  }
  
  # Time averaging step
  averaged_rows <- floor(timesteps_pres / time_window_pres)
  averaged_samples <- vector("list", averaged_rows)
  for (j in 1:averaged_rows) {
    start <- (j - 1) * time_window_pres + 1
    end <- j * time_window_pres
    averaged_samples[[j]] <- as.vector(traitmatrix[start:end, ])
  }
  
  unique_variants <- sort(unique(unlist(averaged_samples)))
  
  # Trait matrix to frequency matrix
  freq_mat <- t(sapply(averaged_samples, function(variants) {
    tab <- table(factor(variants, levels = unique_variants))
    as.numeric(tab) / length(variants)  # relative to N_pres * time_window_pres
  }))
  colnames(freq_mat) <- unique_variants
  
  # Prepare FIT input
  freq_long <- as.data.frame(freq_mat) %>%
    mutate(time = 1:nrow(.)) %>%
    pivot_longer(-time, names_to = "variant", values_to = "freq") %>%
    filter(freq > 0) %>%
    mutate(variant = as.integer(variant))
  
  # Filter to 3+ time points
  freq_long_filtered <- freq_long %>%
    group_by(variant) %>%
    filter(n_distinct(time) >= 3) %>%
    ungroup()
  
  # FIT application
  fit_results <- freq_long_filtered %>%
    group_split(variant) %>%
    map_dfr(~ {
      df <- as.data.frame(.x)
      if (nrow(df) < 3) {
        return(data.frame(variant = df$variant[1], time_points = nrow(df), fit_p = NA, stringsAsFactors = FALSE))
      }
      
      res <- tryCatch(
        fit(time = df$time, v = df$freq),
        error = function(e) list(fit_p = NA)
      )
      
      data.frame(
        variant = df$variant[1],
        time_points = nrow(df),
        fit_p = res$fit_p,
        stringsAsFactors = FALSE
      )
    }) %>%
    mutate(
      sig = ifelse(fit_p > p_value_lvl_pres, "neutral", "selection"),
      sig = ifelse(is.na(fit_p), "NA", sig)
    )
  
  # Store metrics
  total_variants <- nrow(fit_results)
  FPR_pres <- sum(fit_results$sig == "selection") / total_variants
  NDR_pres <- sum(fit_results$sig == "neutral") / total_variants
  
  accuracy_ta_pres[run] <- NDR_pres
}

# Check results
FPR_pres
NDR_pres
accuracy_snapshot[70] # we can check each run individually
overall_accuracy_pres <- mean(accuracy_snapshot, na.rm = TRUE) # mean accuracy across runs
overall_accuracy_pres

# Proportion of the runs have a 95% of detection
over95 <- sum(accuracy_snapshot >= 0.95)
high_accuracy_runs_pres <- over95/n_runs*100
high_accuracy_runs_pres

# Proportion of NA from the simulation
sumNA_pres <- sum(fit_results$sig == "NA")
sumNA_pres
proportionNA_pres <- sumNA_pres/length(fit_results$sig)*100
proportionNA_pres








