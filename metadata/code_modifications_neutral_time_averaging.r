# This is the chunk of code added to emulate time averaging

# Time averaging step -----------------
  averaged_rows <- floor(timesteps / time_window)
  averaged_matrix <- matrix(NA, nrow = averaged_rows, ncol = N) # new matrix
  for (j in 1:averaged_rows) {
    start <- (j - 1) * time_window + 1  # indicate start of window (1, 26, 51...)
    end <- j * time_window # indicate end of window (25, 50, 75...)
    averaged_matrix[j, ] <- apply(traitmatrix[start:end, ], 2, function(x) sample(x, 1))
  }

# And the plot is pretty much the same, adding the alpha to the subtitle and the type of model to the title
ggplot(data.frame(TNR = accuracy_ta), aes(x = TNR)) +
  geom_histogram(binwidth = 0.001, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Neutral Rate (TNR) Across Runs 'Time Averaged' Model", 
       subtitle = "Red line = expected TNR (1 - α)", 
       x = "True Neutral Rate", 
       y = "Frequency",
       caption = paste("Average =", round(mean(accuracy_ta), 3), "|", 
                       "Runs ≥ 95% =", round(high_accuracy_runs, 1), "%", "|",
                       "Number of runs =", n_runs)) +
  theme_minimal()
  
# Stored the plots in an object and then grid both snapshot and time averaged neutral models
plot_neutral_ta <- ggplot(data.frame(TNR = accuracy_ta), aes(x = TNR)) +
  geom_histogram(binwidth = 0.002, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Neutral Rate (TNR) Across Runs 'Time Averaged' Model", 
       subtitle = "Red line = expected TNR (1 - α)", 
       x = "True Neutral Rate", 
       y = "Frequency",
       caption = paste("Average =", round(mean(accuracy_ta), 3), "|", 
                       "Runs ≥ 95% =", round(high_accuracy_runs, 1), "%", "|",
                       "Number of runs =", n_runs)) +
  theme_minimal()

grid.arrange(plot_neutral_snapshot,plot_neutral_ta,ncol=1)

# I modified the caption of the figure from "Average" to "Mean"
plot_neutral_ta <- ggplot(data.frame(TNR = accuracy_ta), aes(x = TNR)) +
  geom_histogram(binwidth = 0.002, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Neutral Rate (TNR) Across Runs 'Time Averaged' Model", 
       subtitle = "Red line = expected TNR (1 - α)", 
       x = "True Neutral Rate", 
       y = "Frequency",
       caption = paste("Mean =", round(mean(accuracy_ta), 3), "|", 
                       "Runs ≥ 95% =", round(high_accuracy_runs, 1), "%", "|",
                       "Number of runs =", n_runs, "|",
                       "% NA =", round(percentageNA, 2), "%")) +
  theme_minimal()
  
# Relevant modification: time averaging went from this:
  averaged_rows <- floor(timesteps / time_window)
  averaged_matrix <- matrix(NA, nrow = averaged_rows, ncol = N) # new matrix
  for (j in 1:averaged_rows) {
    start <- (j - 1) * time_window + 1  # indicate start of window (1, 26, 51...)
    end <- j * time_window # indicate end of window (25, 50, 75...)
    averaged_matrix[j, ] <- apply(traitmatrix[start:end, ], 2, function(x) sample(x, 1))
  }
  
  unique_variants <- sort(unique(as.vector(averaged_matrix))) # store variants
  
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

# To this:
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

# Slightly modified the ggplot function:
plot_neutral_ta <- ggplot(data.frame(NDR = accuracy_ta), aes(x = NDR)) +
  geom_histogram(binwidth = 0.002, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Neutral Detection Rate (NDR) Across Runs 'Time Averaged' Model", 
       subtitle = "Red line = expected NDR (1 - α)", 
       x = "Neutral Detection Rate", 
       y = "Frequency",
       caption = paste("Mean =", round(mean(accuracy_ta), 3), "|", 
                       "Runs ≥ 95% =", round(high_accuracy_runs, 1), "%", "|",
                       "Number of runs =", n_runs, "|",
                       "% NA =", round(proportionNA, 2), "%")) +
  theme_minimal()
 
 # And the "traits" for "variants" in the freq_mat (pipeline, 161-172):
  freq_mat <- t(sapply(averaged_samples, function(variants) {
    tab <- table(factor(variants, levels = unique_variants))
    as.numeric(tab) / length(variants)  # Proportions relative to N * time_window
  }))
  colnames(freq_mat) <- unique_variants
  
# Added a new plot showing p-value distribution
p_value_distribution_ta <- ggplot(data = fit_results, aes(x = fit_p)) +
  geom_histogram(binwidth = 0.025, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 0.075, y = 5, label = "Neutrality", 
           color = "red", size = 3.5, fontface = "bold", hjust = 0) +
  labs(
    title = "P-value Distribution Across Runs 'Time Averaging Model'",
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

p_value_distribution_ta

grid.arrange(p_value_distribution_snapshot, p_value_distribution_ta, ncol = 1)

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
  list(N=100, mu=0.02, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10, time_window = 20),
  list(N=100, mu=0.02, burnin=1000, timesteps=1000, p_value_lvl=0.05, n_runs=10, time_window = 30)
)

all_results_neutral_ta <- map_dfr(params_neutral_ta, ~ {
  sim <- do.call(neutral_ta, args = .x)
  tibble(N  = .x$N,
         mu = .x$mu,
         burnin = .x$burnin,
         timesteps = .x$timesteps,
         p_value_lvl = .x$p_value_lvl,
         n_runs = .x$n_runs,
         time_window = .x$time_window
         mean_accuracy = sim$mean_accuracy,
         high_accuracy_runs = sim$high_accuracy_runs,
         proportionNA = sim$proportionNA,
         mean_p_value = mean(sim$all_pvals, na.rm = TRUE)
  )
})

print(all_results_neutral_ta)

# Multiple innovation rates (20/05):
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


# FOCAL VARIANT PROBLEM SOLUTION 26/05/25
# Record the focal variant as the modal at equilibrium
    pop_counts <- table(pop)
    foc_variant_snap <- as.integer(names(pop_counts)[which.max(pop_counts)])
    
    # Store the focal variant
    results_snap$variant[run] <- foc_variant_snap



# FPR storing (TO ADD TO GITHUB)
accuracy_ta <- numeric(n_runs) # empty vector for accuracy tracking each run
FPR_ta <- numeric(n_runs)
fit_p_count <- vector("list", n_runs) # store p-values

# Total variants = non-NA 
total_variants <- sum(fit_results$sig != "NA")

# FPR calculation
accuracy_ta[run] <- NDR  # track for each run
FPR_ta[run] <- FPR
# Store results
  mean_accuracy <- mean(accuracy_snapshot, na.rm = TRUE) # mean accuracy across runs
  high_accuracy_runs <- sum(accuracy_snapshot >= 0.95) / n_runs * 100
  mean_FPR <- mean(FPR_ta, na.rm = TRUE)

