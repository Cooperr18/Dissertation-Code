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

