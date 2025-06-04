#########################################################
################ NEUTRAL TIME AVERAGING #################
#########################################################

# Install signatselect
install.packages("pak")
pak::pkg_install("benmarwick/signatselect")

# Reading packages
pkgs <- c(
  "signatselect","dplyr","ggplot2",
  "tidyr","gridExtra","purrr",
  "tibble","writexl"
)
lapply(pkgs, library, character.only = TRUE)


# Parameters
N # Number of individuals
mu # innovation rate (numeric, between 0 and 1, inclusive)
burnin # number of initial steps (iterations) discarded
timesteps # actual number of time steps or "generations" after the burn-in
p_value_lvl # Significance level
n_runs # number of test runs
time_window # time window size

accuracy_snapshot <- numeric(n_runs) # empty vector for accuracy tracking each run
fit_p_count <- vector("list", n_runs) # store p-values


# Pipeline as a function ------------------------------
set.seed(1234)

neutral_ta <- function(N, mu, burnin, timesteps, p_value_lvl, n_runs, time_window) {
  
  accuracy_ta <- numeric(n_runs) # empty vector for accuracy tracking each run
  fit_p_count <- vector("list", n_runs) # store p-values
  
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
    
    # Time averaging step
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
    freq_mat <- t(sapply(averaged_samples, function(variants) {
      tab <- table(factor(variants, levels = unique_variants))
      as.numeric(tab) / length(variants)  # Proportions relative to N * time_window
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
    
    # p-value count
    fit_p_count[[run]] <- fit_results$fit_p
    
    # Store metrics
    total_variants <- nrow(fit_results)
    FPR <- sum(fit_results$sig == "selection") / total_variants  # False positives
    NDR <- sum(fit_results$sig == "neutral") / total_variants     # True negatives
    
    if(total_variants == 0) { # if no variants survive NA is returned
      FPR <- NA; NDR <- NA
    }
    
    accuracy_ta[run] <- NDR  # track NDR across runs
  }
  
  # Store results
  mean_accuracy <- mean(accuracy_ta, na.rm = TRUE) # mean accuracy across runs
  high_accuracy_runs <- sum(accuracy_ta >= 0.95) / n_runs * 100
  
  all_pvals <- unlist(fit_p_count)
  
  # Proportion of NA from the simulation
  sumNA <- sum(is.na(all_pvals))
  proportionNA <- sumNA / length(all_pvals) * 100
  if(sumNA == 0) {
    proportionNA <- 0
  }
  
  # Store output
  return(list(
    accuracy_ta = accuracy_ta,
    mean_accuracy = mean_accuracy,
    high_accuracy_runs = high_accuracy_runs,
    sumNA = sumNA,
    proportionNA = proportionNA,
    fit_p_count = fit_p_count,
    N = N,
    mu = mu,
    burnin = burnin,
    timesteps = timesteps,
    p_value_lvl = p_value_lvl,
    n_runs = n_runs,
    time_window = time_window,
    all_pvals = all_pvals
  ))
}


# Run simulation with parameters
n_ta_sim <- neutral_ta(N = 100, mu = 0.02, burnin = 1000,
                             timesteps = 1000, p_value_lvl = 0.05,
                             n_runs = 10, time_window = 20)

# Check output
n_ta_sim

# Store output across runs in a table
results_table_neutral_ta <- tibble(
  N  = n_ta_sim$N,
  mu = n_ta_sim$mu,
  burnin = n_ta_sim$burnin,
  timesteps = n_ta_sim$timesteps,
  p_value_lvl = n_ta_sim$p_value_lvl,
  n_runs = n_ta_sim$n_runs,
  time_window = n_ta_sim$time_window,
  mean_accuracy = n_ta_sim$mean_accuracy,
  high_accuracy_runs = n_ta_sim$high_accuracy_runs,
  proportionNA = n_ta_sim$proportionNA,
  mean_p_value = mean(n_ta_sim$all_pvals, na.rm = TRUE)
)
results_table_neutral_ta

# Run many parameter-sets and stack the results
params_neutral_ta <- list(
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100, time_window = 20),
  list(N=100, mu=0.025, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100, time_window = 20),
  list(N=100, mu=0.05, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100, time_window = 20),
  list(N=100, mu=0.075, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100, time_window = 20),
  list(N=100, mu=0.1, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100, time_window = 20),
  list(N=100, mu=0.125, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100, time_window = 20),
  list(N=100, mu=0.15, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100, time_window = 20),
  list(N=100, mu=0.175, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100, time_window = 20),
  list(N=100, mu=0.2, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100, time_window = 20)
)

all_results_neutral_ta <- map_dfr(params_neutral_ta, ~ {
  sim <- do.call(neutral_ta, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         p_value_lvl = .x$p_value_lvl,
         n_runs = .x$n_runs,
         time_window = .x$time_window,
         mean_accuracy = round(sim$mean_accuracy, 3),
         proportionNA = round(sim$proportionNA, 2),
         mean_p_value = round(mean(sim$all_pvals, na.rm = TRUE), 3)
  )
})

print(all_results_neutral_ta)

write_xlsx(all_results_neutral_snapshot, "neutral_ta_results.xlsx")


# PLOTS ----

# Plot distribution of NDR, marking the 95% threshold
plot_neutral_ta <- function(n_ta_sim, binwidth = 0.005) {
  ggplot(data = data.frame(NDR = n_ta_sim$accuracy_ta), aes(x = NDR)) +
    geom_histogram(binwidth = binwidth, fill = "skyblue", color = "black") +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "Neutral Detection Rate (NDR) Across Runs 'Time Averaged' Model", 
         subtitle = "Red line = expected NDR (1 - α)", 
         x = "Neutral Detection Rate", 
         y = "Frequency",
         caption = paste("Average =", round(mean(n_ta_sim$accuracy_ta), 3), "|", 
                         "Runs ≥ 95% =", round(n_ta_sim$high_accuracy_runs, 1), "%", "|",
                         "Number of runs =", n_ta_sim$n_runs, "|",
                         "% NA =", round(n_ta_sim$proportionNA, 2), "%")) +
    theme_minimal()
}

plot_neutral_ta(n_ta_sim)

# P-value distribution
p_value_distr_ta <- function(n_ta_sim, binwidth = 0.025) {
  ggplot(data = data.frame(p_value = n_ta_sim$all_pvals), aes(x = p_value)) +
    geom_histogram(binwidth = binwidth, fill = "skyblue", color = "black", na.rm = TRUE) +
    geom_vline(xintercept = n_ta_sim$p_value_lvl, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "Distribution of FIT P-values Across Runs Time Averaging Model", 
         subtitle = paste0("Red line = α threshold (", n_ta_sim$p_value_lvl, ")"), 
         x = "p-value", 
         y = "Frequency",
         caption = paste("Average NDR =", round(mean(n_ta_sim$accuracy_ta, na.rm = TRUE), 3), "|",
                         "Mean p-values =", round(mean(n_ta_sim$all_pvals, na.rm = TRUE), 3), "|",
                         "Runs =", n_ta_sim$n_runs, "|",
                         "% NA =", round(n_ta_sim$proportionNA, 2),"%"
         )
    ) +
    theme_minimal()
}

p_value_distr_ta(n_ta_sim)

grid.arrange(p_value_distribution_ta, p_value_distribution_ta, ncol = 1)




