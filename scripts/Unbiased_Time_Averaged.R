#########################################################
################ NEUTRAL TIME AVERAGING #################
#########################################################

# Reading packages
pkgs <- c(
  "pak","dplyr","ggplot2",
  "tidyr","gridExtra","purrr",
  "tibble","writexl"
)
lapply(pkgs, library, character.only = TRUE)

# Install signatselect
pak::pkg_install("benmarwick/signatselect")
library(signatselect)


# PARAMETERS -----------------------------------------------------------------
# N = Population size
# mu = innovation rate (between 0 and 1, inclusive)
# burnin = number of initial steps (iterations) discarded
# timesteps = actual number of time steps or "generations" after the burn-in
# p_value_lvl = Significance level
# n_runs = number of test runs
# time_window = time window/bins size

set.seed(1234)

# PIPELINE -------------------------------------------------------------------

neutral_ta <- function(N, mu, burnin, timesteps, p_value_lvl, n_runs, time_window) {
  
  accuracy_ta <- numeric(n_runs) # empty vector for accuracy tracking each run
  FPR_ta <- numeric(n_runs)
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
    if (nrow(freq_long_filtered)==0) {
      fit_results <- tibble(
        variant     = integer(),
        time_points = integer(),
        fit_p       = double(),
        sig         = character()
      )
    } else {
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
          sig = case_when(
            is.na(fit_p) ~ "NA",
            fit_p > p_value_lvl ~ "neutral",
            TRUE  ~ "selection"
          )
        )
    }
    
    # p-value count
    fit_p_count[[run]] <- fit_results$fit_p
    
    # Store metrics
    total_variants <- sum(fit_results$sig != "NA")
    FPR <- sum(fit_results$sig == "selection") / total_variants  # False positives
    NDR <- sum(fit_results$sig == "neutral") / total_variants     # True negatives
    
    if(total_variants == 0) { # if no variants survive NA is returned
      FPR <- NA; NDR <- NA
    }
    
    accuracy_ta[run] <- NDR  # track NDR across runs
    FPR_ta[run] <- FPR
  }
  
  # Store results
  # MEAN
  mean_accuracy <- mean(accuracy_ta, na.rm = TRUE) # mean accuracy across runs
  high_accuracy_runs <- sum(accuracy_ta >= 0.95) / n_runs * 100
  mean_FPR <- mean(FPR_ta, na.rm = TRUE)
  
  # SD
  sd_NDR <- sd(accuracy_ta, na.rm = T) # sd across runs
  sd_FPR <- sd(FPR_ta, na.rm = T)
  
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
    time_window = time_window,
    
    # OUTPUT
    accuracy_ta = accuracy_ta,
    mean_accuracy = mean_accuracy,
    mean_FPR = mean_FPR,
    sd_NDR = sd_NDR,
    sd_FPR = sd_FPR,
    high_accuracy_runs = high_accuracy_runs,
    sumNA = sumNA,
    proportionNA = proportionNA,
    fit_p_count = fit_p_count,
    all_pvals = all_pvals
    )
  )
}


# BASELINE SIMULATION
n_ta_sim <- neutral_ta(N = 100, mu = 0.02, burnin = 1000,
                             timesteps = 1000, p_value_lvl = 0.05,
                             n_runs = 10, time_window = 20)

# Store output across runs in a table
results_table_neutral_ta <- tibble(
  N  = n_ta_sim$N,
  "µ" = n_ta_sim$mu,
  "B" = n_ta_sim$burnin,
  "Time steps" = n_ta_sim$timesteps,
  "α" = n_ta_sim$p_value_lvl,
  NDR = n_ta_sim$mean_accuracy,
  FPR = n_ta_sim$mean_FPR,
  ">95% runs" = n_ta_sim$high_accuracy_runs,
  "%NA" = n_ta_sim$proportionNA,
  mean_p_value = mean(n_ta_sim$all_pvals, na.rm = TRUE),
  "Runs" = n_ta_sim$n_runs
)
results_table_neutral_ta

