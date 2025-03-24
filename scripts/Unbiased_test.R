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
non_significant_count <- 0

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
  
  # Compute freq as ratio of the count of each variant to the total count at each time step
  freq_long <- freq_long %>%
    group_by(time) %>%
    mutate(total_count = sum(freq)) %>%
    mutate(freq = freq / total_count) %>%
    ungroup()
  
  # Filter out variants that appear in fewer than 3 time points
  freq_long_filtered <- freq_long %>%
    group_by(variant) %>%
    filter(n() >= 3) %>% # fit() requires more than 3 points
    ungroup()
  
  print(head(freq_long_filtered))
  
  fit_safely <- function(time, v) {
    tryCatch({
      fit_result <- fit(time = time, v = v)
      return(list(result = fit_result))  # Return the result as a list
    }, error = function(e) {
      return(list(result = NULL))  # Return NULL in case of an error
    })
  }
  
  df_fit_test_results <- freq_long_filtered %>%
    group_by(variant) %>%
    nest() %>%
    mutate(fit_test = map(data, ~fit_safely(time = .x$time, v = .x$freq))) %>%
    mutate(fit_p = map(fit_test, ~.x$result$p.value)) %>%
    unnest(fit_p) %>%
    mutate(sig = ifelse(fit_p <= 0.05, "selection", "neutral"))
  
  # Checking the results
  head(df_fit_test_results)
  summary(freq_long)
  summary(freq_long_filtered)

  # Run the fit() function on the entire dataset (no filtering)
  fit_result <- tryCatch({
    fit(time = freq_long_filtered$time, v = freq_long_filtered$freq)
  }, error = function(e) {
    return(NULL)  # Handle errors in the fit() function (e.g., if it fails)
  })

  # Check if the p-value is non-significant (p > 0.05)
  if (!is.null(fit_result)) {
    print(fit_result)  # Print the result for debugging
    if (fit_result$p.value > p_value_lvl) {
      non_significant_count <- non_significant_count + 1  # Count non-significant result
    }
  }
}

non_significant_proportion  <- non_significant_count / n_runs
print(paste("Proportion of non-significant results: ", non_significant_proportion))

print(non_significant_count)
print(fit_result)

# --- PLOTS ---

# One single plot
ggplot(freq_long_filtered, aes(x = time, y = freq, group = variant, color = as.factor(variant))) +
  geom_line(alpha = 0.8) +
  labs(title = "Frequency Trajectories of All Variants", x = "Time", y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "none")

# Facet
ggplot(freq_long_filtered, aes(x = time, y = freq, group = variant)) +
  geom_line() +
  facet_wrap(~variant) +
  theme_minimal(base_size = 8) +
  labs(x = "Time", y = "Frequency") +
  theme(legend.position = "none")

                   
                   
                   
                   