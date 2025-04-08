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