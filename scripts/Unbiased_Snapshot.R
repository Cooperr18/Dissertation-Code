#########################################################
################## NEUTRAL SNAPSHOT #####################
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

set.seed(1234)

# PIPELINE -------------------------------------------------------------------

neutral_snapshot <- function(N, mu, burnin, timesteps, p_value_lvl, n_runs) {
  
  accuracy_snapshot <- numeric(n_runs) # empty vector for accuracy tracking each run
  FPR_snapshot <- numeric(n_runs)
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
    total_variants <- sum(fit_results$sig != "NA")
    FPR <- sum(fit_results$sig == "selection") / total_variants  # False positives
    NDR <- sum(fit_results$sig == "neutral") / total_variants     # True negatives
    
    if(total_variants == 0) { # if no variants survive NA is returned
      FPR <- NA; NDR <- NA
    }
    
    accuracy_snapshot[run] <- NDR  # track for each run
    FPR_snapshot[run] <- FPR
  }
  
  # Store results
  # MEAN
  mean_accuracy <- mean(accuracy_snapshot, na.rm = TRUE) # mean accuracy across runs
  high_accuracy_runs <- sum(accuracy_snapshot >= 0.95) / n_runs * 100
  mean_FPR <- mean(FPR_snapshot, na.rm = TRUE)
  
  # SD
  sd_NDR <- sd(accuracy_snapshot, na.rm = T) # sd across runs
  sd_FPR <- sd(FPR_snapshot, na.rm = T)
  
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
    
    # OUTPUTS
    accuracy_snapshot = accuracy_snapshot,
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
n_snap_sim <- neutral_snapshot(N = 100, mu = 0.01, burnin = 1000,
                              timesteps = 1000, p_value_lvl = 0.05, n_runs = 100)

# Store output across runs in a table
results_table_neutral_snapshot <- tibble(
  N  = n_snap_sim$N,
  "µ" = n_snap_sim$mu,
  "B" = n_snap_sim$burnin,
  "Time steps" = n_snap_sim$timesteps,
  "α" = n_snap_sim$p_value_lvl,
  NDR = n_snap_sim$mean_accuracy,
  FPR = n_snap_sim$mean_FPR,
  ">95% runs" = n_snap_sim$high_accuracy_runs,
  "%NA" = n_snap_sim$proportionNA,
  mean_p_value = mean(n_snap_sim$all_pvals, na.rm = TRUE),
  "Runs" = n_snap_sim$n_runs
)
results_table_neutral_snapshot
