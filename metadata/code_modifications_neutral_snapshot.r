# Previous syntax to compute the number of neutral variants detected per run -------
# Observed vs total neutral detection
# count neutral cases per run
neutral_counts_per_run_snapshot[run] <- sum(fit_results$sig == "neutral", na.rm = TRUE)

# count total variants tested and the expected proportion
total_variants_tested <- nrow(fit_results)
expected_neutral_count <- round(0.95 * total_variants_tested) # Expected proportion based on threshold

# observed neutral variants in one run
actual_neutral_count <- sum(fit_results$sig == "neutral", na.rm = TRUE)

# match rate
neutral_match_rate <- actual_neutral_count / expected_neutral_count

# store results per run in the empty object
accuracy_snapshot[run] <- neutral_match_rate

# New Syntax ------------------------------
# Store metrics:
  total_variants <- nrow(fit_results)
  FPR <- sum(fit_results$sig == "selection") / total_variants  # False positives
  TNR <- sum(fit_results$sig == "neutral") / total_variants     # True negatives
  # or
  TNR <- 1 - FPR

  
  accuracy_snapshot[run] <- TNR  # Store true negatives across runs
  
# Also deleted ----------------------------
total_neutral_detections <- sum(neutral_counts_per_run_snapshot)
overall_accuracy <- mean(accuracy_snapshot, na.rm = TRUE) # mean accuracy across runs
# Summary
cat("Total neutral detections:", total_neutral_detections, "\n")
cat("Mean accuracy of FIT test across runs:", overall_accuracy, "\n")
# Both objects are not useful anymore given the current syntax ---------
  
# New plot of TNR distribution ---------------------------
# Plot distribution of TNR, marking the 95% threshold
ggplot(data.frame(TNR = accuracy_snapshot), aes(x = TNR)) +
  geom_histogram(binwidth = 0.009, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Neutral Rate (TNR) Across Runs", 
       subtitle = "Red line = expected TNR (1 - p_value_lvl)", 
       x = "True Neutral Rate", 
       y = "Frequency",
       caption = paste("Average =", round(mean(accuracy_snapshot), 3), "|", 
                       "Runs ≥ 95% =", round(high_accuracy_runs, 1), "%", "|",
                       "Number of runs =", n_runs)) +
  theme_minimal()
  
# Stored the plots in an object and then grid both snapshot and time averaged neutral models
plot_neutral_snapshot <- ggplot(data.frame(TNR = accuracy_snapshot), aes(x = TNR)) +
  geom_histogram(binwidth = 0.005, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Neutral Rate (TNR) Across Runs 'Snapshot' Model", 
       subtitle = "Red line = expected TNR (1 - α)", 
       x = "True Neutral Rate", 
       y = "Frequency",
       caption = paste("Average =", round(mean(accuracy_snapshot), 3), "|", 
                       "Runs ≥ 95% =", round(high_accuracy_runs, 1), "%", "|",
                       "Number of runs =", n_runs)) +
  theme_minimal()

grid.arrange(plot_neutral_snapshot,plot_neutral_ta,ncol=1)

# Slight modification to the plot, from "Average" to "Mean"
plot_neutral_snapshot <- ggplot(data.frame(TNR = accuracy_snapshot), aes(x = TNR)) +
  geom_histogram(binwidth = 0.005, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Neutral Rate (TNR) Across Runs 'Snapshot' Model", 
       subtitle = "Red line = expected TNR (1 - α)", 
       x = "True Neutral Rate", 
       y = "Frequency",
       caption = paste("Mean =", round(mean(accuracy_snapshot), 3), "|", 
                       "Runs ≥ 95% =", round(high_accuracy_runs, 1), "%", "|",
                       "Number of runs =", n_runs)) +
  theme_minimal()
  
# Slight modification to the plot, added the reported % of NA
plot_neutral_snapshot <- ggplot(data.frame(TNR = accuracy_snapshot), aes(x = TNR)) +
  geom_histogram(binwidth = 0.005, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Neutral Rate (TNR) Across Runs 'Snapshot' Model", 
       subtitle = "Red line = expected TNR (1 - α)", 
       x = "True Neutral Rate", 
       y = "Frequency",
       caption = paste("Average =", round(mean(accuracy_snapshot), 3), "|", 
                       "Runs ≥ 95% =", round(high_accuracy_runs, 1), "%", "|",
                       "Number of runs =", n_runs, "|",
                       "% NA =", round(percentageNA, 2), "%")) +
  theme_minimal()

# Slightly modified the ggplot function, from TNR to NDR:
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

# Added a new plot showing p-value distribution:
p_value_distribution_snapshot <- ggplot(data = fit_results, aes(x = fit_p)) +
  geom_histogram(binwidth = 0.025, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 0.075, y = 5, label = "Neutrality", 
           color = "red", size = 3.5, fontface = "bold", hjust = 0) +
  labs(
    title = "P-value Distribution Across Runs 'Snapshot Model'",
    subtitle = paste("µ =", mu, "| N =", N, "| timesteps =", timesteps, "| burnin =", burnin),
    x = "P-value",
    y = "Counts",
    caption = paste(
      "Mean =", round(mean(fit_results$fit_p, na.rm = TRUE), 3), "|",
      "Number of runs =", n_runs, "|",
      "% NA =", round(proportionNA, 2), "%"
    )
  ) +
  coord_cartesian(clip = "off") +  # allows text to overflow if needed
  theme_minimal()

p_value_distribution_snapshot



# Pipeline into function:
neutral_snapshot <- function(N, mu, burnin, timesteps, p_value_lvl, n_runs) {

  neutral_counts_per_run_snapshot <- numeric(n_runs) # empty vector for counting neutral variants
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
    
    # Store metrics for this run
    total_variants <- nrow(fit_results)
    FPR <- sum(fit_results$sig == "selection") / total_variants  # False positives
    NDR <- sum(fit_results$sig == "neutral") / total_variants     # True negatives
    
    accuracy_snapshot[run] <- NDR  # track for each run
  }
  
  # Store results
  mean_accuracy <- mean(accuracy_snapshot, na.rm = TRUE) # mean accuracy across runs
  
  # Proportion runs have 95% detection
  over95 <- sum(accuracy_snapshot >= 0.95)
  high_accuracy_runs <- over95/n_runs*100
  
  # Proportion NA across runs
  sumNA <- sum(fit_results$sig == "NA")
  proportionNA <- sumNA/length(fit_results$sig)*100
  
  # p-value count
  fit_p_count[[run]] <- fit_results$fit_p
  
  # Store output
  return(list(
    accuracy_snapshot = accuracy_snapshot,
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
    n_runs = n_runs
  ))
}

# Run simulation with parameters
n_snap_sim <- neutral_snapshot(N = 100, mu = 0.02, burnin = 1000,
                              timesteps = 1000, p_value_lvl = 0.05, n_runs = 100)

# Functions for both plots:
# NDR
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
  p_vals <- unlist(n_snap_sim$fit_p_count)
  ggplot(data = data.frame(p_value = p_vals), aes(x = p_value)) +
    geom_histogram(binwidth = binwidth, fill = "skyblue", color = "black", na.rm = TRUE) +
    geom_vline(xintercept = n_snap_sim$p_value_lvl, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "DIstribution of FIT P-values Across Runs", 
         subtitle = paste0("Red line = α threshold (", n_snap_sim$p_value_lvl, ")"), 
        x = "p-value", 
        y = "Frequency",
        caption = paste("Total p-values =", length(p_vals), "|",
                        "Runs =", n_snap_sim$n_runs, "|",
                        "% NA =", round(mean(is.na(p_vals)) * 100, 2),"%"
                        )
        ) +
    theme_minimal()
}

p_value_distr_snap(n_snap_sim)

# Additions 17/05/25

fit_p_count <- vector("list", n_runs) # store p-values

fit_p_count[[run]] <- fit_results$fit_p  # store all the p-values from this run


if(total_variants == 0) { # if no variants survive NA is returned
    FPR <- NA; NDR <- NA
  }


all_pvals   <- unlist(fit_p_count) # list of p-values

# Proportion of NA from the simulation
sumNA <- sum(is.na(all_pvals))
sumNA
proportionNA <- sumNA / length(all_pvals) * 100
proportionNA
if(sumNA == 0) {
    proportionNA <- 0
  }

# P-values across runs
mean_p_value <- mean(all_pvals, na.rm = T)
mean_p_value

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

# Run many parameter-sets and stack the results
params_neutral_snapshot <- list(
  list(N=100, mu=0.02, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10),
  list(N=200, mu=0.02, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10)
)

all_results_neutral_snapshot <- map_dfr(params_neutral_snapshot, ~ {
  sim <- do.call(neutral_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         p_value_lvl = .x$p_value_lvl,
         n_runs = .x$n_runs,
         mean_accuracy = sim$mean_accuracy,
         high_accuracy_runs = sim$high_accuracy_runs,
         proportionNA = sim$proportionNA,
         mean_p_value = mean(sim$all_pvals, na.rm = TRUE)
  )
})

print(all_results_neutral_snapshot)

# Multiple innovation rates:
params_neutral_snapshot <- list(
  list(N=100, mu=0.01, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10),
  list(N=100, mu=0.025, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10),
  list(N=100, mu=0.05, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10),
  list(N=100, mu=0.075, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10),
  list(N=100, mu=0.1, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10),
  list(N=100, mu=0.125, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10),
  list(N=100, mu=0.15, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10),
  list(N=100, mu=0.175, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10),
  list(N=100, mu=0.2, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10)
)

all_results_neutral_snapshot <- map_dfr(params_neutral_snapshot, ~ {
  sim <- do.call(neutral_snapshot, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         p_value_lvl = .x$p_value_lvl,
         n_runs = .x$n_runs,
         mean_accuracy = round(sim$mean_accuracy, 3),
         proportionNA = round(sim$proportionNA, 2),
         mean_p_value = round(mean(sim$all_pvals, na.rm = TRUE), 3)
  )
})
print(all_results_neutral_snapshot)

# Export table
write_xlsx(all_results_neutral_snapshot, "neutral_snapshot_results.xlsx")
