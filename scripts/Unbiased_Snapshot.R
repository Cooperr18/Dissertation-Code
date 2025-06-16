#########################################################
################## NEUTRAL SNAPSHOT #####################
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

# PARAMETERS -----------------------------------------------------------------
# N = Population size
# mu = innovation rate (between 0 and 1, inclusive)
# burnin = number of initial steps (iterations) discarded
# timesteps = actual number of time steps or "generations" after the burn-in
# p_value_lvl = Significance level
# n_runs = number of test runs

set.seed(1234)

# PIPELINE -------------------------------------------------------------------

neutral_snapshot <- function(N, mu, burnin, timesteps, p_value_lvl, n_runs) {

  accuracy_snapshot <- numeric(n_runs) # empty vector for accuracy tracking each run 
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
    
    unique_variants <- sort(unique(as.vector(traitmatrix)))
    
    # Trait matrix to frequency matrix
    freq_mat <- t(apply(traitmatrix, 1, function(row) {
      tab <- table(factor(row, levels = unique_variants))
      as.numeric(tab)/N # Convert to frequencies
    }))
    colnames(freq_mat) <- unique_variants # give names
    
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
    
    # Store metrics for this run
    total_variants <- nrow(fit_results)
    FPR <- sum(fit_results$sig == "selection") / total_variants  # False positives
    NDR <- sum(fit_results$sig == "neutral") / total_variants     # True negatives
    
    if(total_variants == 0) { # if no variants survive NA is returned
      FPR <- NA; NDR <- NA
    }
    
    accuracy_snapshot[run] <- NDR  # track for each run
  }
  
  # Store results
  mean_accuracy <- mean(accuracy_snapshot, na.rm = TRUE) # mean accuracy across runs
  high_accuracy_runs <- sum(accuracy_snapshot >= 0.95) / n_runs * 100
  
  all_pvals <- unlist(fit_p_count)
  
  # Proportion of NA from the simulation
  sumNA <- sum(is.na(all_pvals))
  proportionNA <- sumNA / length(all_pvals) * 100
  if(sumNA == 0) {
    proportionNA <- 0
  }
  
  # Store output
  return(list(
    
    # PARAMETERS
    N = N,
    mu = mu,
    burnin = burnin,
    timesteps = timesteps,
    p_value_lvl = p_value_lvl,
    n_runs = n_runs,
    all_pvals = all_pvals,
    
    # OUTPUTS
    accuracy_snapshot = accuracy_snapshot,
    mean_accuracy = mean_accuracy,
    high_accuracy_runs = high_accuracy_runs,
    sumNA = sumNA,
    proportionNA = proportionNA,
    fit_p_count = fit_p_count
    )
  )
}

# Run simulation with one parameter combination
n_snap_sim <- neutral_snapshot(N = 100, mu = 0.02, burnin = 1000,
                              timesteps = 1000, p_value_lvl = 0.05, n_runs = 10)

# Check output
n_snap_sim

# Store output across runs in a table
results_table_neutral_snapshot <- tibble(
  N  = n_snap_sim$N,
  mu = n_snap_sim$mu,
  burnin = n_snap_sim$burnin,
  timesteps = n_snap_sim$timesteps,
  p_value_lvl = n_snap_sim$p_value_lvl,
  n_runs = n_snap_sim$n_runs,
  mean_accuracy = n_snap_sim$mean_accuracy,
  high_accuracy_runs = n_snap_sim$high_accuracy_runs,
  proportionNA = n_snap_sim$proportionNA,
  mean_p_value = mean(n_snap_sim$all_pvals, na.rm = TRUE)
)
results_table_neutral_snapshot


# Run many parameter-sets and stack the results ------------------------------------

# Innovation rate ------------------------------------------------------------------
n_snap_mu_params <- list(
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.025, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.05, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.075, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.1, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.125, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.15, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.175, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.2, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100)
)

# Population size ------------------------------------------------------------------
n_snap_N_params <- list(
  list(N=10, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=50, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=150, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=200, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=250, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=300, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=350, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=400, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100)
)

# Time series ----------------------------------------------------------------------
n_snap_time_params <- list(
  list(N=100, mu=0.01, burnin=100, timesteps=100, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, burnin=200, timesteps=200, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, burnin=500, timesteps=500, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, burnin=750, timesteps=750, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, burnin=1500, timesteps=1500, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, burnin=2000, timesteps=2000, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, burnin=2500, timesteps=2500, p_value_lvl=0.05, n_runs=100),
  list(N=100, mu=0.01, burnin=3000, timesteps=3000, p_value_lvl=0.05, n_runs=100)
)


# Run and store --------------------------------------------------------------------

# INNOVATION RATE (MU)
n_snap_mu_results <- map_dfr(n_snap_mu_params, ~ {
  sim <- do.call(neutral_snaphot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "Time window size" = .x$time_window,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})

# POPULATION SIZE (N)
n_snap_N_results <- map_dfr(n_snap_N_params, ~ {
  sim <- do.call(neutral_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "Time window size" = .x$time_window,
         "α" = .x$p_value_lvl,
         "NDR" = round(sim$mean_accuracy, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Runs" = .x$n_runs
  )
})

# TIME SERIES
n_snap_time_results <- map_dfr(n_snap_time_params, ~ {
  sim <- do.call(neutral_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         "α" = .x$p_value_lvl,
         "Runs" = .x$n_runs,
         "NDR" = round(sim$mean_accuracy, 3),
         "%NA" = round(sim$proportionNA, 2),
         "Mean p-value" = round(mean(sim$all_pvals, na.rm = TRUE), 3)
  )
})

# PRINT RESULTS
print(n_snap_mu_results)
print(n_snap_N_results)
print(n_snap_time_results)

# EXPORT TO SPREADSHEET
write_xlsx(n_snap_mu_results, "tables/n_snap_output/n_snap_mu_params.xlsx") # mu

write_xlsx(n_snap_N_results, "tables/n_snap_output/n_snap_N_params.xlsx") # N

write_xlsx(n_snap_time_params, "tables/n_snap_output/n_snap_time_params.xlsx") # time series



# PLOTS ------------------------------------------------------------------------

# Plot distribution of NDR, marking the 95% threshold
plot_neutral_snapshot <- function(n_snap_sim, binwidth = 0.005) {
  ggplot(data = data.frame(NDR = n_snap_sim$accuracy_snapshot), aes(x = NDR)) +
    geom_histogram(binwidth = binwidth, fill = "skyblue", color = "black") +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "Neutral Detection Rate (NDR) Across Runs 'Snapshot' Model", 
         subtitle = "Red line = expected NDR (1 - α)", 
         x = "Neutral Detection Rate", 
         y = "Frequency",
         caption = paste("Average =", round(mean(n_snap_sim$accuracy_snapshot), 3), "|", 
                         "Runs ≥ 95% =", round(n_snap_sim$high_accuracy_runs, 1), "%", "|",
                         "Number of runs =", n_snap_sim$n_runs, "|",
                         "% NA =", round(n_snap_sim$proportionNA, 2), "%")) +
    theme_minimal()
}

plot_neutral_snapshot(n_snap_sim)

# P-value distribution
p_value_distr_snap <- function(n_snap_sim, binwidth = 0.025) {
  ggplot(data = data.frame(p_value = n_snap_sim$all_pvals), aes(x = p_value)) +
    geom_histogram(binwidth = binwidth, fill = "skyblue", color = "black", na.rm = TRUE) +
    geom_vline(xintercept = n_snap_sim$p_value_lvl, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "Dstribution of FIT P-values Across Runs Snapshot Model", 
         subtitle = paste0("Red line = α threshold (", n_snap_sim$p_value_lvl, ")"), 
        x = "p-value", 
        y = "Frequency",
        caption = paste("Average NDR =", round(mean(n_snap_sim$accuracy_snapshot, na.rm = TRUE), 3), "|",
                        "Mean p-values =", round(mean(n_snap_sim$all_pvals, na.rm = TRUE), 3), "|",
                        "Runs =", n_snap_sim$n_runs, "|",
                        "% NA =", round(n_snap_sim$proportionNA, 2),"%"
                        )
        ) +
    theme_minimal()
}

p_value_distr_snap(n_snap_sim)

# Visualization
freq_long_sig <- freq_long %>%
  left_join(fit_results %>% select(variant, sig), by = "variant")

# One single plot
ggplot(freq_long_filtered, aes(x = time, y = freq, group = variant)) +
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



                   
                   
                   
                   